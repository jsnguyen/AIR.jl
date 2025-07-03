using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LsqFit
using Optim
using ForwardDiff

using AIR
import AIR.crop 

function align_to_template(frame::AstroImage, template::AstroImage; σ::Real=5.0, fillval=NaN)
    # 1) cross‐correlate
    cc = imfilter(frame.data, centered(template.data))

    # 2) find integer peak
    _, ci = findmax(cc)
    r, c = Tuple(ci)

    # 3) sub‐pixel Gaussian fit: p = [A, x0, y0, σ]
    p0  = [cc[ci], Float64(c), Float64(r), σ]
    fit = gaussian2d_fixedwidth_fit(cc, p0, σ)
    x0, y0 = fit[2], fit[3]

    # 4) compute 1‐based center of cc
    cy = (size(cc,1)) / 2
    cx = (size(cc,2)) / 2

    # 5) offset = (Δrow, Δcol)
    offset = (y0 - cy, x0 - cx)
    @info "offset" offset=offset

    # 6) warp the raw data by this sub‐pixel shift
    warped = warp(frame.data,
                  Translation(offset...),
                  axes(frame.data),
                  fill=fillval)

    # 7) return with original header
    return AstroImage(warped, frame.header)
end

"""
    measure_background(frame::AstroImage; mask_radius=50, edge_buffer=20)

Measure the background level in a frame while masking out the PSF.
Returns the median background level from an annular region.
"""
function measure_background(frame::AstroImage; mask_radius=50)
    data = frame.data
    rows, cols = size(data)
    
    # Find the center of the PSF (brightest pixel)
    _, center_idx = findmax(data)
    cy, cx = Tuple(center_idx)
    
    # Create mask to exclude PSF and edges
    mask = trues(size(data))
    
    # Mask out the PSF (circular region around center)
    for i in 1:rows, j in 1:cols
        r = sqrt((i - cy)^2 + (j - cx)^2)
        if r < mask_radius
            mask[i, j] = false
        end
    end
    
    # Extract background pixels
    background_pixels = data[mask]
    
    if length(background_pixels) == 0
        @warn "No background pixels found, returning 0"
        return 0.0
    end
    
    # Return median background level
    return median(background_pixels)
end

"""
    subtract_psf_with_shift(target, reference; initial_shift=(0.0, 0.0), search_radius=5.0, initial_scale=1.0, initial_offset=0.0)

Subtract reference PSF from target PSF after optimizing sub-pixel shift, scaling, and offset using autodiff.
Returns the residual image and optimal parameters (shift, scale, offset).
"""
function subtract_psf_with_shift(target, reference; initial_shift=(0.0, 0.0), search_radius=5.0, initial_scale=1.0, initial_offset=0.0)
    
    # Clean input data - handle NaN/Inf values
    target_clean = replace(target, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    reference_clean = replace(reference, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    
    @info "Input data stats:"
    @info "  Target: $(size(target_clean)), finite elements: $(count(isfinite, target_clean))/$(length(target_clean))"
    @info "  Reference: $(size(reference_clean)), finite elements: $(count(isfinite, reference_clean))/$(length(reference_clean))"
    
    # Objective function: minimize residual after shifting, scaling, and offsetting reference
    function objective(params)
        dy, dx, scale, offset = params[1], params[2], params[3], params[4]
        
        # Skip if shift is too large (enforce hard bounds)
        if abs(dy) > search_radius || abs(dx) > search_radius
            return 1e10
        end
        
        # Skip if scale is unreasonable (between 0.1 and 10.0)
        if scale < 0.1 || scale > 10.0
            return 1e10
        end
        
        try
            # Shift the reference PSF
            shifted_ref = warp(reference_clean, Translation(dy, dx), axes(reference_clean), fill=0.0)
            
            # Apply scaling and offset: scaled_ref = scale * shifted_ref + offset
            scaled_ref = scale .* shifted_ref .+ offset
            
            # Calculate residual (sum of squared differences)
            residual = target_clean .- scaled_ref
            
            # Handle NaN/Inf values
            residual_clean = replace(residual, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
            
            # Calculate sum of squared residuals
            ssr = sum(residual_clean.^2)
            
            # Return large value if result is not finite
            if !isfinite(ssr)
                @debug "Non-finite SSR at params ($dy, $dx, $scale, $offset): $ssr"
                return 1e10
            end
            
            return ssr
            
        catch e
            @debug "Error in objective at params ($dy, $dx, $scale, $offset): $e"
            return 1e10
        end
    end
    
    # Test the objective function at a few points for debugging
    @info "Testing objective function:"
    @info "  At (0,0,1,0): $(objective([0.0, 0.0, 1.0, 0.0]))"
    @info "  At (1,0,1,0): $(objective([1.0, 0.0, 1.0, 0.0]))"
    @info "  At (0,1,1,0): $(objective([0.0, 1.0, 1.0, 0.0]))"
    @info "  At (-1,0,1,0): $(objective([-1.0, 0.0, 1.0, 0.0]))"
    
    # Use bounded LBFGS optimization with explicit bounds
    @info "Starting bounded LBFGS optimization with initial parameters: shift=$(initial_shift), scale=$(initial_scale), offset=$(initial_offset)"
    
    # Set bounds for the optimization [dy, dx, scale, offset]
    lower_bounds = [-search_radius, -search_radius, 0.1, -1000.0]
    upper_bounds = [search_radius, search_radius, 10.0, 1000.0]
    
    result = optimize(objective, 
                     lower_bounds, upper_bounds,
                     [initial_shift[1], initial_shift[2], initial_scale, initial_offset],
                     Fminbox(LBFGS()),
                     Optim.Options(iterations = 200,
                                  g_tol = 1e-6,
                                  show_trace = false))
    
    optimal_params = Optim.minimizer(result)
    final_score = Optim.minimum(result)
    
    @info "Bounded LBFGS optimization completed"
    @info "Final score: $(final_score)"
    @info "Converged: $(Optim.converged(result))"
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

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
    sequences = load_sequences(sequence_obslog)

    as209_4 = load(joinpath(sequences_folder, "as209_4_aligned_frames.fits"), :)
    hbc650_3 = load(joinpath(sequences_folder, "hbc650_3_aligned_frames.fits"), :)

    cube = cat([f.data for f in as209_4]...; dims=3)
    as209_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

    cube = cat([f.data for f in hbc650_3]...; dims=3)
    hbc650_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

    # Perform PSF subtraction with sub-pixel shifting, scaling, and offset
    @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."
    residual, optimal_params = subtract_psf_with_shift(as209_median, hbc650_median, search_radius=20.0)
    
    @info "Optimal shift: ($(round(optimal_params.shift[1], digits=3)), $(round(optimal_params.shift[2], digits=3)))"
    @info "Optimal scale: $(round(optimal_params.scale, digits=3))"
    @info "Optimal offset: $(round(optimal_params.offset, digits=3))"
    @info "Residual statistics: min=$(round(minimum(residual), digits=3)), max=$(round(maximum(residual), digits=3)), std=$(round(std(residual), digits=3))"
    
    # Save results
    save(joinpath(sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequences_folder, "hbc650_median.fits"), hbc650_median)
    save(joinpath(sequences_folder, "psf_subtraction_residual.fits"), residual)
    
    # Also save the optimally shifted, scaled, and offset reference for comparison
    shifted_hbc650 = warp(hbc650_median, Translation(optimal_params.shift...), axes(hbc650_median), fill=0.0)
    shifted_hbc650_clean = replace(shifted_hbc650, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    scaled_shifted_hbc650 = optimal_params.scale .* shifted_hbc650_clean .+ optimal_params.offset
    save(joinpath(sequences_folder, "hbc650_median_shifted.fits"), shifted_hbc650_clean)
    save(joinpath(sequences_folder, "hbc650_median_scaled_shifted.fits"), scaled_shifted_hbc650)



end