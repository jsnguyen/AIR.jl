
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


function fit_generic_kernel(data, initial_guess, kernel)

    function loss(params)
        rows, cols = size(data)
        xs = (1:cols)'
        ys = 1:rows

        residual = data .- kernel.(xs, ys, params...)
        return sum(residual .^ 2)
    end

    res = optimize(loss, initial_guess, LBFGS())
    return Optim.minimizer(res)

end

function fit_2d_gaussian(data, initial_guess; fixed_sigma=nothing, fixed_offset=nothing)

    kernel = if fixed_sigma !== nothing && fixed_offset !== nothing
        if length(initial_guess) != 3
            error("Initial guess must have 3 parameters when both fixed_sigma and fixed_offset are provided")
        end
        (xs, ys, amp, x0, y0) -> gaussian_2d(xs, ys, amp, x0, y0, fixed_sigma, fixed_sigma, fixed_offset)
    elseif fixed_sigma !== nothing
        if length(initial_guess) != 4
            error("Initial guess must have 4 parameters when fixed_sigma is provided")
        end
        (xs, ys, amp, x0, y0, offset) -> gaussian_2d(xs, ys, amp, x0, y0, fixed_sigma, fixed_sigma, offset)
    elseif fixed_offset !== nothing
        if length(initial_guess) != 5
            error("Initial guess must have 5 parameters when fixed_offset is provided")
        end
        (xs, ys, amp, x0, y0, σx, σy) -> gaussian_2d(xs, ys, amp, x0, y0, σx, σy, fixed_offset)
    else
        gaussian_2d
    end

    return fit_generic_kernel(data, initial_guess, kernel)

end

function gaussian_2d(x, y, A, x0, y0, σx, σy, offset)
    return A * exp(-((x - x0)^2 / (2 * σx^2) + (y - y0)^2 / (2 * σy^2))) + offset
end

function fit_and_crop(data, crop_size, initial_guess; fixed_sigma=nothing)

    if fixed_sigma !== nothing
        fit_params = fit_2d_gaussian(data, initial_guess, fixed_sigma=fixed_sigma)
    else
        fit_params = fit_2d_gaussian(data, initial_guess)
    end

    final_cx = fit_params[2] 
    final_cy = fit_params[3]

    data, oy, ox = subpixel_crop(data, crop_size, (final_cy, final_cx))

    return data, final_cx, final_cy, oy, ox

end

function cross_correlation_center(image, template, sigma)
    cc = imfilter(image, centered(template))

    _, max_idx = findmax(image)
    initial_cy, initial_cx = Float64.(Tuple(max_idx))
    initial_guess = [cc[max_idx], initial_cx, initial_cy, median(cc)]
    fit_params = fit_2d_gaussian(cc, initial_guess; fixed_sigma=sigma)

    cx = fit_params[2] 
    cy = fit_params[3]

    return cy, cx

end

function cross_correlate_align(image, template, sigma)

    cy, cx = cross_correlation_center(image, template, sigma)

    offset_x = cx - size(image, 2) / 2
    offset_y = cy - size(image, 1) / 2

    warped = warp(image, Translation(offset_y, offset_x), axes(image))

    return warped

end