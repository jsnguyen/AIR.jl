function gaussian2d_fit(data::Matrix{Float64}, initial_guess::Vector{Float64})
    # Define the 2D Gaussian model
    model(coords, p) = p[1] .* exp.(-((coords[:,1] .- p[2] .- 0.5).^2 ./ (2 * p[4]^2) + (coords[:, 2] .- p[3] .- 0.5).^2 ./ (2 * p[5]^2))) .+ p[6]

    # Create the coordinate grid
    rows, cols = size(data)

    x = repeat(1.0:cols, inner=rows) # x-coordinates repeated for each row
    y = repeat(1.0:rows, outer=cols) # y-coordinates repeated for each column
    coords = hcat(x, y)  # shape is (2, rows*cols)

    # Perform the fit
    fit = curve_fit(model, coords, vec(data), initial_guess)

    return fit.param
end

function gaussian2d_fixedwidth_fit(data::Matrix{Float64}, initial_guess::Vector{Float64}, width::Float64)
    # Define the 2D Gaussian model
    # 1-based indexing here means we have to subtract 0.5 from the x and y coordinates to center the gaussian on the pixel grid
    model(coords, p) = p[1] .* exp.(-((coords[:,1] .- p[2] .- 0.5).^2 ./ (2 * width^2) + (coords[:, 2] .- p[3] .- 0.5).^2 ./ (2 * width^2))) .+ p[4]

    # Create the coordinate grid
    rows, cols = size(data)

    x = repeat(1.0:cols, inner=rows) # x-coordinates repeated for each row
    y = repeat(1.0:rows, outer=cols) # y-coordinates repeated for each column
    coords = hcat(x, y)  # shape is (2, rows*cols)

    # Perform the fit
    fit = curve_fit(model, coords, vec(data), initial_guess)

    return fit.param
end

"""
    fit_gaussian_center_lstsq(data; sigma=5.0)

Fit a 2D Gaussian using least squares to find the PSF center with sub-pixel precision.
Returns a tuple: (center_y, center_x, amplitude, background, fit_quality).

# Arguments
- `data`: 2D array containing the PSF image
- `sigma`: Standard deviation (width) of the Gaussian model (default: 5.0)

# Returns
- `center_y`: Y-coordinate of the fitted center (row index)
- `center_x`: X-coordinate of the fitted center (column index)  
- `amplitude`: Peak amplitude of the fitted Gaussian
- `background`: Fitted background level
- `fit_quality`: Normalized residual metric (lower is better)

The function uses a hybrid approach:
1. Initial guess from brightest pixel location
2. Nonlinear optimization (Nelder-Mead) for center position (cy, cx)
3. Linear least squares for amplitude and background at each trial center

Model: I(i,j) = amplitude * exp(-0.5 * ((i-cy)² + (j-cx)²) / σ²) + background

# Example
```julia
center_y, center_x, amplitude, background, quality = fit_gaussian_center_lstsq(psf_data, sigma=3.0)
```
"""
function fit_gaussian_center_lstsq(data; sigma=5.0)
    rows, cols = size(data)
    
    # Initial guess: brightest pixel location
    _, max_idx = findmax(data)
    initial_cy, initial_cx = Float64.(Tuple(max_idx))
    
    @info "Starting Gaussian center fit with initial guess: ($(round(initial_cy, digits=2)), $(round(initial_cx, digits=2)))"
    
    # Variables to store best fit results
    best_amplitude = 0.0
    best_background = 0.0
    best_residual = Inf
    
    # Objective function: finds optimal center position
    function objective(center_params)
        cy, cx = center_params[1], center_params[2]
        
        # Ensure center is within image bounds
        if cy < 1 || cy > rows || cx < 1 || cx > cols
            return 1e10
        end
        
        # Build linear system for least squares: A * [amplitude, background] = b
        n_pixels = rows * cols
        A = zeros(n_pixels, 2)  # Design matrix
        b = zeros(n_pixels)     # Observed intensities
        
        pixel_idx = 1
        for i in 1:rows, j in 1:cols
            # Gaussian term at pixel (i,j) for center (cy,cx)
            r_squared = (i - cy)^2 + (j - cx)^2
            gaussian_term = exp(-0.5 * r_squared / sigma^2)
            
            A[pixel_idx, 1] = gaussian_term  # Coefficient for amplitude
            A[pixel_idx, 2] = 1.0            # Coefficient for background
            b[pixel_idx] = data[i, j]        # Observed intensity
            pixel_idx += 1
        end
        
        # Solve for amplitude and background using least squares
        try
            fit_params = A \ b
            amplitude, background = fit_params[1], fit_params[2]
            
            # Calculate residual sum of squares
            predicted = A * fit_params
            residuals = b - predicted
            ssr = sum(residuals.^2)
            
            # Penalize unphysical solutions
            if amplitude <= 0
                ssr += 1e10  # Amplitude must be positive
            end
            if background < 0
                ssr += 1e6   # Background should be non-negative
            end
            
            # Update best fit parameters if this is the best so far
            if ssr < best_residual
                best_amplitude = amplitude
                best_background = background
                best_residual = ssr
            end
            
            return ssr
            
        catch e
            @debug "Least squares failed for center ($cy, $cx): $e"
            return 1e10
        end
    end
    
    # Optimize center position using Nelder-Mead simplex method
    @info "Optimizing Gaussian center position..."
    result = optimize(objective, 
                     [initial_cy, initial_cx],
                     NelderMead(),
                     Optim.Options(iterations = 200, 
                                  f_tol = 1e-8))
    
    # Extract optimization results
    optimal_center = Optim.minimizer(result)
    final_residual = Optim.minimum(result)
    converged = Optim.converged(result)
    
    cy_fit, cx_fit = optimal_center[1], optimal_center[2]
    
    # Calculate normalized fit quality metric
    total_signal = sum(abs.(data))
    fit_quality = sqrt(best_residual / length(data)) / (total_signal / length(data))
    
    @info "Gaussian fit completed:"
    @info "  Center: ($(round(cy_fit, digits=3)), $(round(cx_fit, digits=3)))"
    @info "  Amplitude: $(round(best_amplitude, digits=1))"
    @info "  Background: $(round(best_background, digits=1))"
    @info "  Fit quality: $(round(fit_quality, digits=4)) (lower is better)"
    @info "  Converged: $converged"
    
    return cy_fit, cx_fit, best_amplitude, best_background, fit_quality
end


"""
    subtract_psf_with_shift(target, reference; initial_shift=(0.0, 0.0), search_radius=5.0, initial_scale=1.0, initial_offset=0.0)

Subtract reference PSF from target PSF after optimizing sub-pixel shift, scaling, and offset using autodiff.
Returns the residual image and optimal parameters (shift, scale, offset).
"""
function subtract_psf_with_shift(target, reference; initial_shift=(0.0, 0.0), search_radius=5.0, initial_scale=1.0, initial_offset=0.0, mask_radius=30)
    
    target_clean = replace(target, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    reference_clean = replace(reference, NaN => 0.0, Inf => 0.0, -Inf => 0.0)

    @info "Input data stats:"
    @info "  Target: $(size(target_clean)), finite elements: $(count(isfinite, target_clean))/$(length(target_clean))"
    @info "  Reference: $(size(reference_clean)), finite elements: $(count(isfinite, reference_clean))/$(length(reference_clean))"
    
    function objective(params)
        dy, dx, scale, offset = params[1], params[2], params[3], params[4]
        
        if abs(dy) > search_radius || abs(dx) > search_radius
            @warn "Shift out of bounds: dy=$dy, dx=$dx (max search radius: $search_radius)"
            return 1e10
        end

        if scale < 0.1 || scale > 10.0
            @warn "Scale out of bounds: scale=$scale (must be between 0.1 and 10.0)"
            return 1e10
        end

        shifted_ref = warp(reference_clean, Translation(dy, dx), axes(reference_clean), fill=0.0)
        scaled_ref = scale .* shifted_ref .+ offset
        residual = target_clean .- scaled_ref

        if mask_radius != 0.0
            circle_mask = make_circle_mask(size(residual), mask_radius)
            residual[circle_mask] .= 0.0  # Mask out central region
        end


        residual_clean = replace(residual, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
        ssr = sum(residual_clean.^2)/ length(residual_clean)
        
        if !isfinite(ssr)
            @warn "Non-finite SSR at params ($dy, $dx, $scale, $offset)"
            return 1e10
        end
        
        return ssr
            
    end
    
    # Use bounded LBFGS optimization with explicit bounds
    @info "Starting bounded LBFGS optimization with initial parameters: shift=$(initial_shift), scale=$(initial_scale), offset=$(initial_offset)"
    
    # Set bounds for the optimization [dy, dx, scale, offset]
    lower_bounds = [-search_radius, -search_radius, 0.1, -1000.0]
    upper_bounds = [search_radius, search_radius, 20.0, 1000.0]
    
    result = optimize(objective, 
                     lower_bounds, upper_bounds,
                     [initial_shift[1], initial_shift[2], initial_scale, initial_offset],
                     Fminbox(LBFGS()),
                     Optim.Options(iterations = 100,
                                  g_tol = 1e-6,
                                  show_trace = false))
    
    optimal_params = Optim.minimizer(result)
    final_score = Optim.minimum(result)
    @info "Final score" final_score=final_score
    
    @info "Bounded LBFGS optimization completed"
    @info "Iterations: $(Optim.iterations(result))"
    @info "Optimal parameters: dy=$(round(optimal_params[1], digits=3)), dx=$(round(optimal_params[2], digits=3)), scale=$(round(optimal_params[3], digits=3)), offset=$(round(optimal_params[4], digits=3))"
    
    # Ensure the optimization actually improved the result
    initial_params = [initial_shift[1], initial_shift[2], initial_scale, initial_offset]
    initial_score = objective(initial_params)
    if final_score >= initial_score
        @warn "Optimization did not improve result (initial: $initial_score, final: $final_score)"
        @warn "Using initial parameters as result"
        optimal_params = initial_params
    end
    
    # Apply optimal shift, scale, and offset, then subtract
    shifted_reference = warp(reference_clean, Translation(optimal_params[1], optimal_params[2]), axes(reference_clean), fill=0.0)
    
    # Clean the shifted reference to handle any NaN/Inf values from interpolation
    shifted_reference_clean = replace(shifted_reference, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    
    # Apply scaling and offset
    scaled_reference = optimal_params[3] .* shifted_reference_clean .+ optimal_params[4]
    
    residual = target_clean .- scaled_reference
    
    # Clean the final residual
    residual_clean = replace(residual, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    
    # Return residual and all optimal parameters as a named tuple
    optimal_result = (
        shift = (optimal_params[1], optimal_params[2]),
        scale = optimal_params[3],
        offset = optimal_params[4]
    )
    
    return residual_clean, optimal_result
end

"""
    fit_gaussian_center_variable_sigma(data; initial_sigma=5.0, min_sigma=1.0, max_sigma=20.0)

Fit a 2D Gaussian using least squares to find the PSF center with sub-pixel precision and variable width.
Returns a tuple: (center_y, center_x, amplitude, background, sigma, fit_quality).

# Arguments
- `data`: 2D array containing the PSF image
- `initial_sigma`: Initial guess for the Gaussian width (default: 5.0)
- `min_sigma`: Minimum allowed sigma value (default: 1.0)
- `max_sigma`: Maximum allowed sigma value (default: 20.0)

# Returns
- `center_y`: Y-coordinate of the fitted center (row index)
- `center_x`: X-coordinate of the fitted center (column index)  
- `amplitude`: Peak amplitude of the fitted Gaussian
- `background`: Fitted background level
- `sigma`: Fitted standard deviation (width) of the Gaussian
- `fit_quality`: Normalized residual metric (lower is better)

The function uses a hybrid approach:
1. Initial guess from brightest pixel location and provided sigma
2. Nonlinear optimization (Nelder-Mead) for center position (cy, cx) AND sigma
3. Linear least squares for amplitude and background at each trial center/sigma

Model: I(i,j) = amplitude * exp(-0.5 * ((i-cy)² + (j-cx)²) / σ²) + background

# Example
```julia
center_y, center_x, amplitude, background, sigma, quality = fit_gaussian_center_variable_sigma(psf_data, initial_sigma=3.0)
```
"""
function fit_gaussian_center_variable_sigma(data; initial_sigma=5.0, min_sigma=1.0, max_sigma=20.0)
    rows, cols = size(data)
    
    # Initial guess: brightest pixel location
    _, max_idx = findmax(data)
    initial_cy, initial_cx = Float64.(Tuple(max_idx))
    
    @info "Starting variable-sigma Gaussian fit with initial guess: center=($(round(initial_cy, digits=2)), $(round(initial_cx, digits=2))), sigma=$(round(initial_sigma, digits=2))"
    
    # Variables to store best fit results
    best_amplitude = 0.0
    best_background = 0.0
    best_sigma = initial_sigma
    best_residual = Inf
    
    # Objective function: finds optimal center position AND sigma
    function objective(params)
        cy, cx, sigma = params[1], params[2], params[3]
        
        # Ensure center is within image bounds
        if cy < 1 || cy > rows || cx < 1 || cx > cols
            return 1e10
        end
        
        # Ensure sigma is within reasonable bounds
        if sigma < min_sigma || sigma > max_sigma
            return 1e10
        end
        
        # Build linear system for least squares: A * [amplitude, background] = b
        n_pixels = rows * cols
        A = zeros(n_pixels, 2)  # Design matrix
        b = zeros(n_pixels)     # Observed intensities
        
        pixel_idx = 1
        for i in 1:rows, j in 1:cols
            # Gaussian term at pixel (i,j) for center (cy,cx) and width sigma
            r_squared = (i - cy)^2 + (j - cx)^2
            gaussian_term = exp(-0.5 * r_squared / sigma^2)
            
            A[pixel_idx, 1] = gaussian_term  # Coefficient for amplitude
            A[pixel_idx, 2] = 1.0            # Coefficient for background
            b[pixel_idx] = data[i, j]        # Observed intensity
            pixel_idx += 1
        end
        
        # Solve for amplitude and background using least squares
        try
            fit_params = A \ b
            amplitude, background = fit_params[1], fit_params[2]
            
            # Calculate residual sum of squares
            predicted = A * fit_params
            residuals = b - predicted
            ssr = sum(residuals.^2)
            
            # Penalize unphysical solutions
            if amplitude <= 0
                ssr += 1e10  # Amplitude must be positive
            end
            if background < 0
                ssr += 1e6   # Background should be non-negative
            end
            
            # Update best fit parameters if this is the best so far
            if ssr < best_residual
                best_amplitude = amplitude
                best_background = background
                best_sigma = sigma
                best_residual = ssr
            end
            
            return ssr
            
        catch e
            @debug "Least squares failed for center ($cy, $cx), sigma=$sigma: $e"
            return 1e10
        end
    end
    
    # Optimize center position AND sigma using Nelder-Mead simplex method
    @info "Optimizing Gaussian center position and width..."
    result = optimize(objective, 
                     [initial_cy, initial_cx, initial_sigma],
                     NelderMead(),
                     Optim.Options(iterations = 300,  # More iterations for 3D optimization
                                  f_tol = 1e-8))
    
    # Extract optimization results
    optimal_params = Optim.minimizer(result)
    final_residual = Optim.minimum(result)
    converged = Optim.converged(result)
    
    cy_fit, cx_fit, sigma_fit = optimal_params[1], optimal_params[2], optimal_params[3]
    
    # Calculate normalized fit quality metric
    total_signal = sum(abs.(data))
    fit_quality = sqrt(best_residual / length(data)) / (total_signal / length(data))
    
    @info "Variable-sigma Gaussian fit completed:"
    @info "  Center: ($(round(cy_fit, digits=3)), $(round(cx_fit, digits=3)))"
    @info "  Sigma: $(round(sigma_fit, digits=3))"
    @info "  Amplitude: $(round(best_amplitude, digits=1))"
    @info "  Background: $(round(best_background, digits=1))"
    @info "  Fit quality: $(round(fit_quality, digits=4)) (lower is better)"
    @info "  Converged: $converged"
    
    return cy_fit, cx_fit, best_amplitude, best_background, sigma_fit, fit_quality
end