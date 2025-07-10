using Printf
using AstroImages

using AIR
import AIR.crop 

@autolog begin

    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
    sequences = load_sequences(sequence_obslog)

    as209 = load(joinpath(sequences_folder, "as209_median.fits"))
    template_psf = load(joinpath(sequences_folder, "as209_1_template_psf.fits"))

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

    save(joinpath(sequences_folder, "as209_cropped.fits"), cropped_target) 
    save(joinpath(sequences_folder, "as209_cropped_psf.fits"), cropped_psf) 

end