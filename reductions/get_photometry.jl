using AstroImages

using AIR
import AIR.crop 

autolog("$(@__FILE__).log") do


    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
    sequences = load_sequences(sequence_obslog)

    as209 = load(joinpath(sequences_folder, "as209_median.fits"))
    template_psf = load(joinpath(sequences_folder, "as209_1_template_psf.fits"))
    
    boxsize = 30
    coarse_location = (128,161)
    cropped = crop(as209, (boxsize, boxsize), center=coarse_location)
    cropped_psf = crop(template_psf, (boxsize-5, boxsize-5))

    cx, cy, _, _, _, _ = fit_gaussian_center_variable_sigma(cropped)

    # Calculate precise distance and position angle from image center
    # Image center coordinates (1-indexed)
    image_center_y = size(as209, 1) / 2 + 0.5  # Row center
    image_center_x = size(as209, 2) / 2 + 0.5  # Column center
    
    # Transform fitted center from cropped coordinates to full image coordinates
    # The cropped image was centered at coarse_location, so we need to add the offset
    # from the crop center to the fitted center position
    
    # Crop center in full image coordinates
    crop_center_y, crop_center_x = coarse_location[1], coarse_location[2]
    
    # Size of the cropped region
    crop_half_y = boxsize ÷ 2
    crop_half_x = boxsize ÷ 2
    
    # Cropped image center coordinates (within the crop)
    cropped_center_y = crop_half_y + 0.5
    cropped_center_x = crop_half_x + 0.5
    
    # Offset of fitted center from crop center (in cropped coordinates)
    offset_from_crop_center_y = cy - cropped_center_y
    offset_from_crop_center_x = cx - cropped_center_x
    
    # Transform to full image coordinates
    object_y = crop_center_y + offset_from_crop_center_y
    object_x = crop_center_x + offset_from_crop_center_x
    
    # Calculate offset from full image center (in pixels)
    # Note: dy positive = south (increasing row), dx positive = east (increasing column)
    # For astronomical convention: dy positive = north, so we need to negate
    dy = object_y - image_center_y  # North-South offset (north is positiv)
    dx = object_x - image_center_x     # East-West offset (east is positive)

    
    # Calculate distance from center (in pixels)
    distance_pixels = sqrt(dy^2 + dx^2)
    
    # Calculate position angle (PA) in degrees
    # PA = 0° = North, 90° = East, 180° = South, 270° = West
    # atan2(dx, dy) gives angle from north (y-axis) towards east (x-axis)
    pa_radians = atan(dx, dy)  # atan2 equivalent: atan(y, x) = atan2(y, x)
    pa_degrees = rad2deg(pa_radians)
    
    # Ensure PA is in range [0, 360)
    if pa_degrees < 0
        pa_degrees += 360.0
    end
    
    @info "Coordinate transformation details:"
    @info "  Full image size: $(size(as209))"
    @info "  Full image center: ($(round(image_center_y, digits=1)), $(round(image_center_x, digits=1)))"
    @info "  Crop center (coarse location): ($crop_center_y, $crop_center_x)"
    @info "  Crop size: $(boxsize) × $(boxsize)"
    @info "  Fitted center in crop: ($(round(cy, digits=3)), $(round(cx, digits=3)))"
    @info "  Offset from crop center: ($(round(offset_from_crop_center_y, digits=3)), $(round(offset_from_crop_center_x, digits=3)))"
    @info ""
    @info "Object position relative to full image center:"
    @info "  Object location (full image): ($(round(object_y, digits=3)), $(round(object_x, digits=3)))"
    @info "  Offset (dy, dx): ($(round(dy, digits=2)), $(round(dx, digits=2))) pixels"
    @info "  Distance from center: $(round(distance_pixels*NIRC2_plate_scale, digits=3)) arcsec"
    @info "  Position angle: $(round(pa_degrees, digits=3))° (0°=N, 90°=E, 180°=S, 270°=W)"

    save(joinpath(sequences_folder, "as209_cropped.fits"), cropped) 
    save(joinpath(sequences_folder, "as209_cropped_psf.fits"), cropped_psf) 



end