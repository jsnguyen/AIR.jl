using Printf
using AstroImages

using AIR
import AIR.crop 

function get_photometry_epoch_1()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"))
    template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_template_psf.fits"))

    boxsize = 30
    template_psf_size = 20
    coarse_location = (128,161)

    cropped, cr_y , cr_x = crop(as209.data, (boxsize, boxsize), center=coarse_location)
    cropped_psf,_ ,_ = crop(template_psf.data, (template_psf_size, template_psf_size))

    initial_guess = [20.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, -2]
    fit_params = fit_2d_gaussian(cropped, initial_guess)

    @info "Fit Params:"
    @info @sprintf("  amp    [e-/s] -> %8.3f", fit_params[1])
    @info @sprintf("  x,y      [px] -> %8.3f , %8.3f", fit_params[2], fit_params[3])
    @info @sprintf("  σx,σy    [px] -> %8.3f , %8.3f", fit_params[4], fit_params[5])
    @info @sprintf("  offset [e-/s] -> %8.3f", fit_params[6])

    object_y = fit_params[3] + cr_y 
    object_x = fit_params[2] + cr_x

    image_center_y = size(as209, 1) / 2 + 0.5
    image_center_x = size(as209, 2) / 2 + 0.5
    dy = object_y - image_center_y  # North-South offset (north is positive)
    dx = object_x - image_center_x  # East-West offset (east is positive)
    distance_pixels = sqrt(dy^2 + dx^2)

    # Calculate position angle (PA) in degrees
    pa = rad2deg(atan(dx, dy))
    if pa < 0
        pa += 360.0
    end

    @info "Object position:"
    @info @sprintf("  Abs Pos (y,x)    [px] -> %8.3f, %8.3f", object_y, object_x)
    @info @sprintf("  Rel Pos (δy, δx) [px] -> %8.3f, %8.3f", dy, dx)

    @info "Object params:"
    @info @sprintf("  Sep [\"] -> %8.3f", distance_pixels*NIRC2_plate_scale)
    @info @sprintf("  PA  [°] -> %8.3f", pa)

    cropped_target, oy, ox = subpixel_crop(as209.data, (boxsize, boxsize), (object_y, object_x))

    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped.fits"), cropped_target) 
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_psf.fits"), cropped_psf) 
end

function get_photometry_epoch_2()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-08-02_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"))
    template_psf = load("data/2002-06-16/sequences/as209_template_psf_cored.fits")

    boxsize = 30
    template_psf_size = 20
    coarse_location = (82, 113)

    cropped, cr_y , cr_x = crop(as209.data, (boxsize, boxsize), center=coarse_location)
    cropped_psf,_ ,_ = crop(template_psf.data, (template_psf_size, template_psf_size))

    initial_guess = [20.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, -2]
    fit_params = fit_2d_gaussian(cropped, initial_guess)

    @info "Fit Params:"
    @info @sprintf("  amp    [e-/s] -> %8.3f", fit_params[1])
    @info @sprintf("  x,y      [px] -> %8.3f , %8.3f", fit_params[2], fit_params[3])
    @info @sprintf("  σx,σy    [px] -> %8.3f , %8.3f", fit_params[4], fit_params[5])
    @info @sprintf("  offset [e-/s] -> %8.3f", fit_params[6])

    object_y = fit_params[3] + cr_y 
    object_x = fit_params[2] + cr_x

    image_center_y = size(as209, 1) / 2 + 0.5
    image_center_x = size(as209, 2) / 2 + 0.5
    dy = object_y - image_center_y  # North-South offset (north is positive)
    dx = object_x - image_center_x  # East-West offset (east is positive)
    distance_pixels = sqrt(dy^2 + dx^2)

    # Calculate position angle (PA) in degrees
    pa = rad2deg(atan(dx, dy))
    if pa < 0
        pa += 360.0
    end

    @info "Object position:"
    @info @sprintf("  Abs Pos (y,x)    [px] -> %8.3f, %8.3f", object_y, object_x)
    @info @sprintf("  Rel Pos (δy, δx) [px] -> %8.3f, %8.3f", dy, dx)

    @info "Object params:"
    @info @sprintf("  Sep [\"] -> %8.3f", distance_pixels*NIRC2_plate_scale)
    @info @sprintf("  PA  [°] -> %8.3f", pa)

    cropped_target, oy, ox = subpixel_crop(as209.data, (boxsize, boxsize), (object_y, object_x))

    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped.fits"), cropped_target) 
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_psf.fits"), cropped_psf) 
end

function get_photometry_epoch_3()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-08-21_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"))
    template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_template_psf.fits"))

    boxsize = 30
    template_psf_size = 20
    coarse_location = (48, 79)

    cropped, cr_y , cr_x = crop(as209.data, (boxsize, boxsize), center=coarse_location)
    cropped_psf,_ ,_ = crop(template_psf.data, (template_psf_size, template_psf_size))

    initial_guess = [20.0, boxsize/2 + 0.5, boxsize/2 + 0.5, 3.0, 3.0, -2]
    fit_params = fit_2d_gaussian(cropped, initial_guess)

    @info "Fit Params:"
    @info @sprintf("  amp    [e-/s] -> %8.3f", fit_params[1])
    @info @sprintf("  x,y      [px] -> %8.3f , %8.3f", fit_params[2], fit_params[3])
    @info @sprintf("  σx,σy    [px] -> %8.3f , %8.3f", fit_params[4], fit_params[5])
    @info @sprintf("  offset [e-/s] -> %8.3f", fit_params[6])

    object_y = fit_params[3] + cr_y 
    object_x = fit_params[2] + cr_x

    image_center_y = size(as209, 1) / 2 + 0.5
    image_center_x = size(as209, 2) / 2 + 0.5
    dy = object_y - image_center_y  # North-South offset (north is positive)
    dx = object_x - image_center_x  # East-West offset (east is positive)
    distance_pixels = sqrt(dy^2 + dx^2)

    # Calculate position angle (PA) in degrees
    pa = rad2deg(atan(dx, dy))
    if pa < 0
        pa += 360.0
    end

    @info "Object position:"
    @info @sprintf("  Abs Pos (y,x)    [px] -> %8.3f, %8.3f", object_y, object_x)
    @info @sprintf("  Rel Pos (δy, δx) [px] -> %8.3f, %8.3f", dy, dx)

    @info "Object params:"
    @info @sprintf("  Sep [\"] -> %8.3f", distance_pixels*NIRC2_plate_scale)
    @info @sprintf("  PA  [°] -> %8.3f", pa)

    cropped_target, oy, ox = subpixel_crop(as209.data, (boxsize, boxsize), (object_y, object_x))

    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped.fits"), cropped_target) 
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cropped_psf.fits"), cropped_psf) 
end

@autolog begin

    get_photometry_epoch_1()
    get_photometry_epoch_2()
    get_photometry_epoch_3()

end