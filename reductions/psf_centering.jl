using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LsqFit

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

function gaussian2d_fit(data::Matrix{Float64}, initial_guess::Vector{Float64})
    # Define the 2D Gaussian model
    model(coords, p) = p[1] .* exp.(-((coords[:,1] .- p[2]).^2 ./ (2 * p[4]^2) + (coords[:, 2] .- p[3]).^2 ./ (2 * p[5]^2))) .+ p[6]

    # Create the coordinate grid
    rows, cols = size(data)

    x = repeat(1.0:cols, inner=rows) # x-coordinates repeated for each row
    y = repeat(1.0:rows, outer=cols) # y-coordinates repeated for each column
    coords = hcat(x, y)  # shape is (2, rows*cols)

    # Perform the fit
    fit = curve_fit(model, coords, vec(data), initial_guess)

    return fit.param
end

function gaussian2d_fixedwidth_fit(data::Matrix{Float64}, initial_guess::Vector{Float64}, width::Float64)
    # Define the 2D Gaussian model
    model(coords, p) = p[1] .* exp.(-((coords[:,1] .- p[2]).^2 ./ (2 * width^2) + (coords[:, 2] .- p[3]).^2 ./ (2 * width^2))) .+ p[4]

    # Create the coordinate grid
    rows, cols = size(data)

    x = repeat(1.0:cols, inner=rows) # x-coordinates repeated for each row
    y = repeat(1.0:rows, outer=cols) # y-coordinates repeated for each column
    coords = hcat(x, y)  # shape is (2, rows*cols)

    # Perform the fit
    fit = curve_fit(model, coords, vec(data), initial_guess)

    return fit.param
end

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    unsat = load_frames(sequence_obslog, "unsat")
    j_seq = load_frames(sequence_obslog, "j_seq")
    kp_seq = load_frames(sequence_obslog, "kp_seq")

    coarse_size = 120
    fine_size = 100
    cropped_frames = AstroImage[]

    for frame in unsat
        _, coarse_center = findmax(frame.data)
        cropped = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
        push!(cropped_frames, AstroImage(cropped))
        @info "coarse center" coarse_center
    end

    template = cropped_frames[1].data
    #template_offset = calculate_centroid_offset(template)

    # want a really tight gaussian for finding the center
    template_gauss_fit_width = 3.0
    template_fit_params = gaussian2d_fixedwidth_fit(template, [1.0e9, size(template, 2)/2, size(template, 1)/2, 1.0e8], template_gauss_fit_width)
    template_offset = (template_fit_params[3] - size(template, 1)/2, template_fit_params[2] - size(template, 2)/2)

    @info "template params" template_offset=template_offset
    template = warp(template, Translation(template_offset...), axes(template)) |> x -> crop(x, (fine_size, fine_size))
    template_offset = calculate_centroid_offset(template)
    @info "template params" template_offset=template_offset

    n_iter = 3
    for i in 1:n_iter

        ccrs = AstroImage[]
        offsets = NTuple{2, Float64}[]

        @info "Calculating sub-pixel cross-correlation for all frames against the template..."
        @showprogress "Correlating frames..." for i in eachindex(cropped_frames)


            cross_corr = imfilter(cropped_frames[i], centered(template))
            push!(ccrs, cross_corr)

            # Find the integer index of the peak in the cross-correlation map.
            _, ind_peak_ccr = findmax(cross_corr)

            gaussian_width = 5.0
            fit_params = gaussian2d_fixedwidth_fit(cross_corr.data, [1.0e9, Float64(ind_peak_ccr[2])-0.5, Float64(ind_peak_ccr[1])-0.5, 1.0e7], gaussian_width)

            peak_x = fit_params[2]
            peak_y = fit_params[3]
            center_y, center_x = size(cross_corr) ./ 2
            offset = (peak_y - center_y, peak_x - center_x) .+ 1.0

            push!(offsets, offset)

            aligned_frames = AstroImage[]
            for (frame, offset) in zip(cropped_frames, offsets)
                translated_frame = warp(frame, Translation(offset...), axes(frame)) |> x -> crop(x, (fine_size, fine_size))
                push!(aligned_frames, AstroImage(translated_frame, frame.header))
            end

            template = median(framelist_to_cube(aligned_frames), dims=3) |> x -> dropdims(x, dims=3)

            sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
            save(joinpath(sequences_folder, "template.fits"), template)
            save(joinpath(sequences_folder, "reg.fits"), framelist_to_cube(cropped_frames))
            save(joinpath(sequences_folder, "aligned.fits"), framelist_to_cube(aligned_frames))
            save(joinpath(sequences_folder, "ccr.fits"), framelist_to_cube(ccrs))

        end

        @info "Finished cross-correlation. Found sub-pixel offsets for $(length(offsets)) frames." offsets=offsets

    end


end