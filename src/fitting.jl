
function optimal_subtract_target(target::AbstractArray, reference::AbstractArray, initial_guess::Vector{Float64}, search_radius::Float64; scale_bounds=(0.1, 10), offset_bounds=(-100.0, 100.0), inner_mask_radius=nothing, outer_mask_radius=nothing)

    function align_and_subtract(params)
        dy, dx, scale, offset = Float64.(params)

        shifted_ref = warp(reference, Translation(dy, dx), axes(reference), fill=0.0)

        scaled_ref = scale .* shifted_ref .+ offset
        residual = target .- scaled_ref

        remove_nan!(residual)

        return residual
    end

    function loss(params)
        residual = align_and_subtract(params)

        inner_circle_mask = trues(size(residual))
        if inner_mask_radius !== nothing && inner_mask_radius > 0.0
            if inner_mask_radius > 0.0
                inner_circle_mask = .~make_circle_mask(size(residual), inner_mask_radius)
            end
        end

        outer_circle_mask = trues(size(residual))
        if outer_mask_radius !== nothing
            if 0.0 < inner_mask_radius < outer_mask_radius
                outer_circle_mask = make_circle_mask(size(residual), outer_mask_radius)
            end
        end
        
        mask = inner_circle_mask .& outer_circle_mask

        annulus = residual[mask]

        lsq = sum(annulus.^2)/ length(annulus)
        return lsq
    end

    lower_bounds = [-search_radius, -search_radius, scale_bounds[1], offset_bounds[1]]
    upper_bounds = [search_radius, search_radius, scale_bounds[2], offset_bounds[2]]

    res = optimize(loss, lower_bounds, upper_bounds, initial_guess, Fminbox(NelderMead()))
    params = Optim.minimizer(res)
    residual = align_and_subtract(params)

    return residual, params

end

function optimal_subtract_target(target::AstroImage, reference::AstroImage, initial_guess::Vector{Float64}, search_radius::Float64; scale_bounds=(0.1, 10), offset_bounds=(-100.0, 100.0), inner_mask_radius=nothing, outer_mask_radius=nothing)

    residual, params = optimal_subtract_target(target.data, reference.data, initial_guess, search_radius; scale_bounds=scale_bounds, offset_bounds=offset_bounds, inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

    return AstroImage(residual, target.header), params
end

function fit_generic_kernel(data, initial_guess, kernel; lower_bounds=nothing, upper_bounds=nothing, bounded_fit=false)

    function loss(params)
        rows, cols = size(data)
        xs = (1:cols)'
        ys = 1:rows

        residual = data .- kernel.(xs, ys, params...)
        return sum(residual .^ 2)
    end

    if bounded_fit
        if lower_bounds === nothing || upper_bounds === nothing
            error("For bounded fit, both lower_bounds and upper_bounds must be provided")
        end
        res = optimize(loss, lower_bounds, upper_bounds, initial_guess, Fminbox(LBFGS()))
    else
        res = optimize(loss, initial_guess, LBFGS())
    end

    return Optim.minimizer(res)

end

function fit_2d_gaussian(data, initial_guess; fixed_sigma=nothing, fixed_offset=nothing, kwargs...)

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

    return fit_generic_kernel(data, initial_guess, kernel; kwargs...)

end

function gaussian_2d(x, y, A, x0, y0, σx, σy, offset)
    return A * exp(-((x - x0)^2 / (2 * σx^2) + (y - y0)^2 / (2 * σy^2))) + offset
end

function gaussian_2d_rotated(x, y, amp, x0, y0, σx, σy, θ, offset)
    cos_θ = cos(θ)
    sin_θ = sin(θ)
    x_rot = cos_θ * (x - x0) + sin_θ * (y - y0)
    y_rot = -sin_θ * (x - x0) + cos_θ * (y - y0)
    return amp * exp(-((x_rot)^2 / (2 * σx^2) + (y_rot)^2 / (2 * σy^2))) + offset
end

function fit_and_crop(data, crop_size, initial_guess; fixed_sigma=nothing)

    if fixed_sigma !== nothing
        fit_params = fit_2d_gaussian(data, initial_guess; fixed_sigma=fixed_sigma)
    else
        fit_params = fit_2d_gaussian(data, initial_guess)
    end

    final_cx = fit_params[2] 
    final_cy = fit_params[3]

    data, oy, ox = subpixel_crop(data, crop_size, (final_cy, final_cx))

    return data, final_cx, final_cy, oy, ox

end

function cross_correlation_center(image::AbstractArray, template::AbstractArray, sigma)
    cc = imfilter(image, centered(template))

    _, max_idx = findmax(image)
    initial_cy, initial_cx = Float64.(Tuple(max_idx))
    initial_guess = [cc[max_idx], initial_cx, initial_cy, median(cc)]
    fit_params = fit_2d_gaussian(cc, initial_guess; fixed_sigma=sigma)

    cx = fit_params[2] 
    cy = fit_params[3]

    return cy, cx

end

function cross_correlate_align(image::AbstractArray, template::AbstractArray, sigma)

    cy, cx = cross_correlation_center(image, template, sigma)

    offset_x = cx - size(image, 2) / 2
    offset_y = cy - size(image, 1) / 2

    warped = warp(image, Translation(offset_y, offset_x), axes(image))

    return warped

end

function cross_correlate_align(image::AstroImage, template::AstroImage, sigma)

    # warp messes up the AstroImage type i think
    warped = cross_correlate_align(image.data, template.data, sigma)
    return AstroImage(warped, image.header)
end