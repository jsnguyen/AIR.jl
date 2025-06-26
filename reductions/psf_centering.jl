using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations

using AIR

"""
Calculates the center of mass (centroid) of a small image patch.
Returns the subpixel offset (y, x) from the center of the patch.
The patch should have odd dimensions for a well-defined center.
"""
function calculate_centroid_offset(patch::AbstractMatrix)
    # Ensure we don't divide by zero
    total_flux = sum(patch)
    total_flux == 0 && return (0.0, 0.0)

    # Create coordinate grids relative to the center of the patch
    h, w = size(patch)
    y_coords = (1:h) .- (h + 1) / 2
    x_coords = (1:w) .- (w + 1) / 2

    # Calculate the weighted average of coordinates (centroid)
    y_offset = sum(patch .* y_coords) / total_flux
    x_offset = sum(patch .* x_coords') / total_flux

    return (y_offset, x_offset)
end

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    unsat = load_frames(sequence_obslog, "unsat")
    j_seq = load_frames(sequence_obslog, "j_seq")
    kp_seq = load_frames(sequence_obslog, "kp_seq")

    cropped_size = 60
    cropped_frames = []

    for frame in unsat
        _, coarse_center = findmax(frame.data)
        cropped = crop(frame, (cropped_size, cropped_size), center=Tuple(coarse_center))
        push!(cropped_frames, cropped)
        @info "coarse center" coarse_center

    end

    template = cropped_frames[1]

    # We'll store the floating point (y, x) offsets for each frame relative to the template.
    offsets = NTuple{2, Float64}[]
    centroid_box_size = 5 # Must be an odd number

    ccrs = AstroImage[]

    @info "Calculating sub-pixel cross-correlation for all frames against the template..."
    @showprogress "Correlating frames..." for i in eachindex(cropped_frames)

        cross_corr = imfilter(padarray(cropped_frames[i], Pad(:reflect, cropped_size, cropped_size)), template)
        push!(ccrs, cross_corr)

        # Find the integer index of the peak in the cross-correlation map.
        _, peak_index_int = findmax(cross_corr)

        # Create a small cutout around the peak to calculate the centroid for sub-pixel precision.
        half_box = centroid_box_size ÷ 2
        # Define the window, ensuring it doesn't go out of bounds
        y_range = max(1, peak_index_int[1] - half_box):min(size(cross_corr, 1), peak_index_int[1] + half_box)
        x_range = max(1, peak_index_int[2] - half_box):min(size(cross_corr, 2), peak_index_int[2] + half_box)
        
        patch = cross_corr[y_range, x_range]

        # Calculate the sub-pixel offset within the patch
        y_sub, x_sub = calculate_centroid_offset(patch)

        # The final sub-pixel peak location is the integer peak plus the centroid offset
        peak_y = peak_index_int[1] + y_sub
        peak_x = peak_index_int[2] + x_sub

        # The center of the map corresponds to zero shift. The offset is the displacement from the center.
        center_y, center_x = size(cross_corr) ./ 2 .+ 0.5
        offset = (peak_y - 2*center_y, peak_x - 2*center_x)
        push!(offsets, offset)
    end

    @info "Finished cross-correlation. Found sub-pixel offsets for $(length(offsets)) frames." offsets=offsets

    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")

    aligned_frames = AstroImage[]
    for (frame, offset) in zip(cropped_frames, offsets)
        translation = Translation(.-offset...) 
        warped_frames = warp(frame, translation, axes(frame))
        push!(aligned_frames, AstroImage(warped_frames, frame.header))

    end
    save(joinpath(sequences_folder, "unsat.fits"), framelist_to_cube(aligned_frames))
    save(joinpath(sequences_folder, "ccr.fits"), framelist_to_cube(ccrs))

end