using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
import CoordinateTransformations: recenter
using LsqFit
using Optim
using ForwardDiff
using Rotations

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
    rotate_image_center(img::AstroImage, angle_degrees; fillval=0.0)

Rotate an image about its center by the specified angle in degrees.
"""
function rotate_image_center(img::AstroImage, angle_degrees; fillval=0.0)
    angle_rad = deg2rad(angle_degrees)
    
    rows, cols = size(img.data)
    center = (rows + 1) / 2, (cols + 1) / 2
    
    # Create rotation about specified center point
    rotation = recenter(RotMatrix(angle_rad), center)
    
    rotated_data = warp(img.data, rotation, axes(img.data), fill=fillval)
    
    return AstroImage(rotated_data, img.header)
end

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
    sequences = load_sequences(sequence_obslog)

    as209_3 = load(joinpath(sequences_folder, "as209_3_aligned_frames.fits"), :)
    as209_4 = load(joinpath(sequences_folder, "as209_4_aligned_frames.fits"), :)

    #as209 = cat(as209_3, as209_4; dims=3)
    as209 = as209_4

    hbc650 = load(joinpath(sequences_folder, "hbc650_3_aligned_frames.fits"), :)

    for frame in hbc650
        frame ./= frame["ITIME"]  # Normalize by integration time
    end
    cube = framelist_to_cube(hbc650)
    hbc650_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

    inner_mask_radius = 60
    outer_mask_radius = 230
    inner_circle_mask = make_circle_mask(size(as209[1]), inner_mask_radius)
    outer_circle_mask = make_circle_mask(size(as209[1]), outer_mask_radius)

    # Perform PSF subtraction with sub-pixel shifting, scaling, and offset
    derotated = AstroImage[]
    for frame in as209
        frame ./= frame["ITIME"]

        @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."
        residual, optimal_params = subtract_psf_with_shift(frame, hbc650_median; mask_radius=inner_mask_radius)

        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle)
        derot[inner_circle_mask] .= NaN
        derot[.!outer_circle_mask] .= NaN
        push!(derotated, derot)
        
        @info "Optimal shift: ($(round(optimal_params.shift[1], digits=3)), $(round(optimal_params.shift[2], digits=3)))"
        @info "Optimal scale: $(round(optimal_params.scale, digits=3))"
        @info "Optimal offset: $(round(optimal_params.offset, digits=3))"
        @info "Residual statistics: min=$(round(minimum(derot), digits=3)), max=$(round(maximum(derot), digits=3)), std=$(round(std(derot), digits=3))"
    end

    as209_median = framelist_to_cube(derotated) |> x -> median(x, dims=3) |> x -> dropdims(x, dims=3)
    as209_median[inner_circle_mask] .= NaN
    as209_median[.!outer_circle_mask] .= NaN
    
    # Save results
    save(joinpath(sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequences_folder, "as209_psf_subtraction_derotated.fits"), derotated...)

end