using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Optim

using AIR
import AIR.crop 



@autolog begin

    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)
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

    search_radius = 5
    inner_mask_radius = 60
    outer_mask_radius = 230
    inner_circle_mask = make_circle_mask(size(as209[1]), inner_mask_radius)
    outer_circle_mask = make_circle_mask(size(as209[1]), outer_mask_radius)

    # Perform PSF subtraction with sub-pixel shifting, scaling, and offset
    derotated = AstroImage[]
    initial_guess = [0.0, 0.0, 1.0, 0.0]
    for frame in as209
        frame ./= frame["ITIME"]

        @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."
        #residual, optimal_params = subtract_psf_with_shift(frame, hbc650_median; mask_radius=inner_mask_radius)
        residual, optimal_params = optimal_subtract_target(frame, hbc650_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)


        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle)
        derot[inner_circle_mask] .= NaN
        derot[.!outer_circle_mask] .= NaN
        push!(derotated, derot)

        initial_guess = optimal_params  # Use the last optimal parameters as the initial guess for the next frame
        @info "Optimal Params" optimal_params=optimal_params

    end

    as209_median = framelist_to_cube(derotated) |> x -> median(x, dims=3) |> x -> dropdims(x, dims=3)
    as209_median[inner_circle_mask] .= NaN
    as209_median[.!outer_circle_mask] .= NaN

    # Save results
    save(joinpath(sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequences_folder, "as209_psf_subtraction_derotated.fits"), derotated...)

end