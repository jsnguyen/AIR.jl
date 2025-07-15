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

    err_sep = sqrt((x/sep)^2 * sigma_x + (y/sep)^2 * sigma_y)
    err_ang = sqrt((y/sep^2)^2 * sigma_x + (x/sep^2)^2 * sigma_y)

    return sep, ang, err_sep, err_ang

end

function get_astrometry(as209_northup_median, as209_northup_cube, template_psf, median_coarse_location, frame_coarse_locations; initial_guess=nothing, boxsize=14, fixed_sigma=3.0)

    if initial_guess === nothing
        #initial_guess = [30.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, 0.0, -2.0]
        initial_guess = [30.0, boxsize/2 + 0.5, boxsize/2 + 0.5, -2.0]
    end

    median_cropped, cr_y , cr_x = crop(as209_northup_median, (boxsize, boxsize), center=median_coarse_location)
    cropped_psf,_ ,_ = crop(template_psf, (boxsize, boxsize))

    #lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, -2*π, -1000.0]
    #upper_bounds = [100.0, boxsize+1, boxsize+1, 10.0, 10.0, 2*π, 1000.0]

    lower_bounds = [0.0, 0.0, 0.0, -1000.0]
    upper_bounds = [100.0, boxsize+1, boxsize+1, 1000.0]

    jitter_xs, jitter_ys = [], []
    as209_northup_cropped_cube = AstroImage[]
    for (coarse_location,frame) in zip(frame_coarse_locations,as209_northup_cube)
        cropped, cr_y , cr_x = crop(frame, (boxsize, boxsize), center=coarse_location)
        #fit_params = fit_2d_gaussian(cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
        #fit_params = fit_generic_kernel(cropped, initial_guess, gaussian_2d_rotated; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
        fit_params = fit_2d_gaussian(cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true, fixed_sigma=fixed_sigma)

        object_y = fit_params[3] + cr_y 
        object_x = fit_params[2] + cr_x

        cropped_frame, oy, ox = subpixel_crop(frame, (boxsize, boxsize), (object_y, object_x))
        push!(as209_northup_cropped_cube, cropped_frame)

        push!(jitter_xs, object_x)
        push!(jitter_ys, object_y)
    end

    std_xs = std(jitter_xs)
    std_ys = std(jitter_ys)

    #fit_params = fit_2d_gaussian(median_cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
    #fit_params = fit_generic_kernel(median_cropped, initial_guess, gaussian_2d_rotated; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)
    fit_params = fit_2d_gaussian(median_cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true, fixed_sigma=fixed_sigma)

    @info "Jitter STD:"
    @info @sprintf("  x [px] -> %8.6f", std_xs)
    @info @sprintf("  y [px] -> %8.6f", std_ys)

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

    image_center_y = size(as209_northup_median, 1) / 2 + 0.5
    image_center_x = size(as209_northup_median, 2) / 2 + 0.5

    sep, pa, err_sep, err_pa = error_cartesian_to_polar(object_x, object_y, image_center_x, image_center_y, std_xs, std_ys)

    sep = sep * NIRC2_plate_scale  # Convert pixels to arcseconds
    err_sep = err_sep * NIRC2_plate_scale  # Convert pixels to arcseconds

    # Calculate position angle (PA) in degrees
    pa = rad2deg(pa)
    if pa < 0
        pa += 360.0
    end

    @info "Object position:"
    @info @sprintf("  Abs Pos (y,x)    [px] -> %8.3f, %8.3f", object_y, object_x)

    @info "Object params:"
    @info @sprintf("  Sep [\"] -> %8.3f +/- %4.3f", sep, err_sep)
    @info @sprintf("  PA  [°] -> %8.3f +/- %4.3f", pa, err_pa)

    cropped_target, oy, ox = subpixel_crop(as209_northup_median, (boxsize, boxsize), (object_y, object_x))

    return cropped_target, cropped_psf, sep, pa, err_sep, err_pa, as209_northup_cropped_cube

end

@autolog begin

    epochs = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]
    median_coarse_locations = [
        (80, 110),  # 2002-06-16
        (77, 107),  # 2002-08-02
        (77, 111),  # 2002-08-21
        (75, 117)   # 2005-07-27
    ]

    frame_coarse_locations = [((80,111),(81,110), (79,109), (80,109)),
                              ((75,108), (76,109), (77,107),(76,108),(78,107),(78,107)),
                              ((77,111), (77, 111), (77,111), (77,111)),
                              ((76,118), (76,119), (75,118), (75,118), (76,117), (76,117), (75,118), (76,117))]


    astrometry = OrderedDict{String, Any}()
    for (i, date) in enumerate(epochs)
        @info "Processing date: $date"

        sequence_obslog = Obslog("pipeline/obslogs/$(date)_sequences.toml")

        as209_northup_median = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_median.fits"))
        as209_northup_cube = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_cube.fits"),:)

        # this date doesnt have a template psf, so we use the one from epoch 1
        if date == "2002-08-02"
            template_psf = load("data/2002-06-16/sequences/as209_template_psf.fits")
        else
            template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_template_psf.fits"))
        end

        cropped_target, cropped_psf, sep, pa, err_sep, err_pa, as209_northup_cropped_cube = get_astrometry(as209_northup_median, as209_northup_cube, template_psf, median_coarse_locations[i], frame_coarse_locations[i])

        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_target.fits"), cropped_target) 
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_template_psf.fits"), cropped_psf) 
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_cropped_cube.fits"), as209_northup_cropped_cube...) 

        astrometry["epoch_$(date)"] = OrderedDict(
            "date" => date,
            "sep" => sep,
            "pa" => pa,
            "err_sep" => err_sep,
            "err_pa" => err_pa
        )

    end

    write_toml(joinpath("pipeline/as209_astrometry.toml"), astrometry)

end