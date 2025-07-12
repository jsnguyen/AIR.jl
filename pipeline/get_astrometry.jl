using Printf
using AstroImages

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

function get_astrometry(as209, as209_cube, template_psf, coarse_location; initial_guess=nothing, boxsize=30, distance_pixel_cutoff=5.0)

    if initial_guess === nothing
        initial_guess = [30.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, -2]
    end

    median_cropped, cr_y , cr_x = crop(as209.data, (boxsize, boxsize), center=coarse_location)
    cropped_psf,_ ,_ = crop(template_psf.data, (boxsize, boxsize))

    lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, -10.0]
    upper_bounds = [100.0, boxsize+1, boxsize+1, 10.0, 10.0, 10.0]

    jitter_xs, jitter_ys = [], []
    for frame in as209_cube
        cropped, cr_y , cr_x = crop(frame.data, (boxsize, boxsize), center=coarse_location)
        fit_params = fit_2d_gaussian(cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)

        crop_center_y = size(cropped, 1) / 2 + 0.5
        crop_center_x = size(cropped, 2) / 2 + 0.5
        dist = sqrt((fit_params[2]-crop_center_x)^2 + (fit_params[3]-crop_center_y)^2)
        if dist > distance_pixel_cutoff
            @info "Skipping frame with bad fit: $(dist) $(fit_params[2]), $(fit_params[3])"
            continue
        end

        push!(jitter_xs, fit_params[2])
        push!(jitter_ys, fit_params[3])
    end

    std_xs = std(jitter_xs)
    std_ys = std(jitter_ys)

    fit_params = fit_2d_gaussian(median_cropped, initial_guess; lower_bounds=lower_bounds, upper_bounds=upper_bounds, bounded_fit=true)

    @info "Jitter STD:"
    @info @sprintf("  x [px] -> %8.6f", std_xs)
    @info @sprintf("  y [px] -> %8.6f", std_ys)

    @info "Fit Params:"
    @info @sprintf("  amp    [e-/s] -> %8.3f", fit_params[1])
    @info @sprintf("  x,y      [px] -> %8.3f , %8.3f", fit_params[2], fit_params[3])
    @info @sprintf("  σx,σy    [px] -> %8.3f , %8.3f", fit_params[4], fit_params[5])
    @info @sprintf("  offset [e-/s] -> %8.3f", fit_params[6])

    object_y = fit_params[3] + cr_y 
    object_x = fit_params[2] + cr_x

    image_center_y = size(as209, 1) / 2 + 0.5
    image_center_x = size(as209, 2) / 2 + 0.5

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

    cropped_target, oy, ox = subpixel_crop(as209.data, (boxsize, boxsize), (object_y, object_x))

    return cropped_target, cropped_psf, sep, pa

end

@autolog begin

    epochs = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]
    coarse_locations = [
        (130, 159),  # 2002-06-16
        (80, 112),   # 2002-08-02
        (48, 79),    # 2002-08-21
        (52, 90)     # 2005-07-27
    ]

    for (date, coarse_location) in zip(epochs, coarse_locations)
        @info "Processing date: $date"

        sequence_obslog = Obslog("pipeline/obslogs/$(date)_sequences.toml")

        as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"))
        as209_cube = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_cube.fits"),:)

        # this date doesnt have a template psf, so we use the one from epoch 1
        if date == "2002-08-02"
            template_psf = load("data/2002-06-16/sequences/as209_template_psf_cored.fits")
        else
            template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_template_psf.fits"))
        end

        cropped_target, cropped_psf, sep, pa = get_astrometry(as209, as209_cube, template_psf, coarse_location)

        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_target.fits"), cropped_target) 
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_template_psf.fits"), cropped_psf) 

    end

end