using Printf
using Statistics
using AstroImages
using Photometry

using AIR
import AIR.crop 

# obtained from error propagation formula
function error_cartesian_to_polar(x, y, center_x, center_y, sigma_x, sigma_y)
    dx = x - center_x
    dy = y - center_y

    sep = sqrt(dy^2 + dx^2)
    ang = atan(dx,dy)

    err_sep = sqrt((x/sep * sigma_x)^2 + (y/sep * sigma_y)^2)

    # note the swap here in y and x
    err_ang = sqrt((y/sep^2 * sigma_x)^2 + (x/sep^2 * sigma_y)^2)

    return sep, ang, err_sep, err_ang

end

function get_astrometry(northup_median, northup_cube, template_psf, median_coarse_location, frame_coarse_locations; initial_guess=nothing, boxsize=14, fixed_sigma=3.0)

    if initial_guess === nothing
        #initial_guess = [30.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, 0.0, -2.0]
        initial_guess = [30.0, boxsize/2 + 0.5, boxsize/2 + 0.5, -2.0]
    end

    median_cropped, cr_y , cr_x = crop(northup_median, (boxsize, boxsize), center=median_coarse_location)
    cropped_psf,_ ,_ = crop(template_psf, (boxsize, boxsize))

    #lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, -2*π, -1000.0]
    #upper_bounds = [100.0, boxsize+1, boxsize+1, 10.0, 10.0, 2*π, 1000.0]

    lower_bounds = [0.0, 0.0, 0.0, -1000.0]
    upper_bounds = [100.0, boxsize+1, boxsize+1, 1000.0]

    jitter_xs, jitter_ys = [], []
    northup_cropped_cube = AstroImage[]
    for (coarse_location,frame) in zip(frame_coarse_locations,northup_cube)
        cropped, cr_y , cr_x = crop(frame, (boxsize, boxsize), center=coarse_location)
        #fit_params = fit_2d_gaussian(cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
        #fit_params = fit_generic_kernel(cropped, initial_guess, gaussian_2d_rotated; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
        fit_params = fit_2d_gaussian(cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true, fixed_sigma=fixed_sigma)

        object_y = fit_params[3] + cr_y 
        object_x = fit_params[2] + cr_x

        cropped_frame, oy, ox = subpixel_crop(frame, (boxsize, boxsize), (object_y, object_x))
        push!(northup_cropped_cube, cropped_frame)

        push!(jitter_xs, object_x)
        push!(jitter_ys, object_y)
    end

    err_xs = std(jitter_xs)
    err_ys = std(jitter_ys)

    @info "Jitter STD:"
    @info @sprintf("  x [px] -> %8.6f", err_xs)
    @info @sprintf("  y [px] -> %8.6f", err_ys)

    err_xs = err_xs*NIRC2_plate_scale * 1000
    err_ys = err_ys*NIRC2_plate_scale * 1000

    #fit_params = fit_2d_gaussian(median_cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
    #fit_params = fit_generic_kernel(median_cropped, initial_guess, gaussian_2d_rotated; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
    fit_params = fit_2d_gaussian(median_cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true, fixed_sigma=fixed_sigma)

    @info @sprintf("  x [mas] -> %8.6f", err_xs)
    @info @sprintf("  y [mas] -> %8.6f", err_ys)

    @info "Fit Params:"
    @info @sprintf("  amp    [e-/s] -> %8.3f", fit_params[1])
    @info @sprintf("  x,y      [px] -> %8.3f , %8.3f", fit_params[2], fit_params[3])
    @info @sprintf("  offset [e-/s] -> %8.3f", fit_params[4])
    #@info @sprintf("  σx,σy    [px] -> %8.3f , %8.3f", fit_params[4], fit_params[5])
    #@info @sprintf("  offset [e-/s] -> %8.3f", fit_params[6])
    #@info @sprintf("  rot     [rad] -> %8.3f", fit_params[6])
    #@info @sprintf("  offset [e-/s] -> %8.3f", fit_params[7])

    object_y = fit_params[3] + cr_y 
    object_x = fit_params[2] + cr_x

    image_center_y = size(northup_median, 1) / 2 + 0.5
    image_center_x = size(northup_median, 2) / 2 + 0.5

    sep, pa, err_sep, err_pa = error_cartesian_to_polar(object_x, object_y, image_center_x, image_center_y, err_xs, err_ys)

    sep = sep * NIRC2_plate_scale  # Convert pixels to arcseconds
    err_sep = err_sep * NIRC2_plate_scale  # Convert pixels to arcseconds

    # Calculate position angle (PA) in degrees
    pa = rad2deg(pa)
    if pa < 0
        pa += 360.0
    end
    err_pa = rad2deg(err_pa)

    @info "Object position:"
    @info @sprintf("  Abs Pos (y,x)    [px] -> %8.3f, %8.3f", object_y, object_x)

    @info "Object params:"
    @info @sprintf("  Sep [\"] -> %8.3f +/- %4.3f", sep, err_sep)
    @info @sprintf("  PA  [°] -> %8.3f +/- %4.3f", pa, err_pa)

    cropped_target, oy, ox = subpixel_crop(northup_median, (boxsize, boxsize), (object_y, object_x))

    return cropped_target, cropped_psf, sep, pa, err_sep, err_pa, err_xs, err_ys, northup_cropped_cube

end

@stage function epoch_get_astrometry(date, median_coarse_location, frame_coarse_locations, northup_median, northup_cube, template_psf; injected_psf=nothing, save_path=nothing)

    if save_path === nothing
        save_path = joinpath("analyze_AS_209/$(date)_astrometry.toml")
    end

    paths = context["paths"]

    if injected_psf !== nothing
        template_psf = injected_psf
    end

    cropped_target, cropped_psf, sep, pa, err_sep, err_pa, err_xs, err_ys, northup_cropped_cube = get_astrometry(northup_median, northup_cube, template_psf, median_coarse_location, frame_coarse_locations)

    save(joinpath(paths.sequences_folder, "$(date)_cropped_target.fits"), cropped_target) 
    save(joinpath(paths.sequences_folder, "$(date)_cropped_template_psf.fits"), cropped_psf) 
    save(joinpath(paths.sequences_folder, "$(date)_northup_cropped_cube.fits"), northup_cropped_cube...) 

    astrometry = OrderedDict{String, Any}()
    astrometry["epoch_$(date)"] = OrderedDict("date" => date,
                                                "sep" => sep,
                                                "sep_units" => "arcsec",
                                                "pa" => pa,
                                                "pa_units" => "degrees",
                                                "err_sep" => err_sep,
                                                "err_sep_units" => "arcsec",
                                                "err_pa" => err_pa,
                                                "err_pa_units" => "degrees",
                                                "err_xs" => err_xs,
                                                "err_xs_units" => "mas",
                                                "err_ys" => err_ys,
                                                "err_ys_units" => "mas"
    )

    write_toml(save_path, astrometry)

    return northup_cropped_cube, template_psf

end