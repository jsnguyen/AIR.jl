using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Optim
using ADI

using AIR
import AIR.crop 

function rdi_subtract_epoch_1()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209_3 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_3_aligned_frames.fits"), :)
    as209_4 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_4_aligned_frames.fits"), :)

    #as209 = cat(as209_3, as209_4; dims=3)
    as209 = as209_4

    hbc650 = load(joinpath(sequence_obslog.paths.sequences_folder, "hbc650_3_aligned_frames.fits"), :)

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
        residual, optimal_params = optimal_subtract_target(frame, hbc650_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle)
        derot[inner_circle_mask] .= NaN
        derot[.!outer_circle_mask] .= NaN
        push!(derotated, AstroImage(derot,frame.header))

        initial_guess = optimal_params  # Use the last optimal parameters as the initial guess for the next frame
        @info "Optimal Params" optimal_params=optimal_params

    end

    as209_median = framelist_to_cube(derotated) |> x -> median(x, dims=3) |> x -> dropdims(x, dims=3)
    as209_median[inner_circle_mask] .= NaN
    as209_median[.!outer_circle_mask] .= NaN

    # Save results
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cube.fits"), derotated...)
end

function rdi_subtract_epoch_2()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-08-02_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209_2 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_2_aligned_frames.fits"), :)
    as209 = as209_2

    derotated = AstroImage[]
    for frame in as209
        frame ./= frame["ITIME"]
        derot = rotate_image_center(copy(frame), -frame["ROTPPOSN"])
        nan_mask = isnan.(derot.data) .| isinf.(derot.data) # create a mask for NaN and Inf values
        derot[nan_mask] .= 0.0 # Set NaN and Inf values to NaN in the derotated frame
        push!(derotated, AstroImage(derot, frame.header))
    end

    median_frame = median(framelist_to_cube(derotated), dims=3) |> x -> dropdims(x, dims=3)

    search_radius = 5
    inner_mask_radius = 60
    outer_mask_radius = 230

    derotated = AstroImage[]
    initial_guess = [0.0, 0.0, 1.0, 0.0]
    for frame in as209
        @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."

        re_rotated_median = rotate_image_center(median_frame, frame["ROTPPOSN"])
        residual, optimal_params = optimal_subtract_target(frame, re_rotated_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

        angle, _ = calculate_north_angle(frame.header)
        residual = rotate_image_center(residual, -angle+frame["PARANG"]) # need parang so the target doesn't move
        push!(derotated, AstroImage(residual, frame.header))

        initial_guess = optimal_params
        @info "Optimal Params" optimal_params=optimal_params
    end

    subtracted_median_frame = median(framelist_to_cube(derotated), dims=3) |> x -> dropdims(x, dims=3)
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"), subtracted_median_frame)
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cube.fits"), derotated...)
end

function rdi_subtract_epoch_3()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-08-21_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209_4 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_4_aligned_frames.fits"), :)

    as209 = as209_4

    tyc2307 = load(joinpath(sequence_obslog.paths.sequences_folder, "tyc2307_3_aligned_frames.fits"), :)

    for frame in tyc2307
        frame ./= frame["ITIME"]  # Normalize by integration time
    end
    cube = framelist_to_cube(tyc2307)
    tyc2307_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

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
        residual, optimal_params = optimal_subtract_target(frame, tyc2307_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle+residual["PARANG"])
        derot[inner_circle_mask] .= NaN
        derot[.!outer_circle_mask] .= NaN
        push!(derotated, AstroImage(derot,frame.header))

        initial_guess = optimal_params  # Use the last optimal parameters as the initial guess for the next frame
        @info "Optimal Params" optimal_params=optimal_params

    end

    as209_median = framelist_to_cube(derotated) |> x -> median(x, dims=3) |> x -> dropdims(x, dims=3)
    as209_median[inner_circle_mask] .= NaN
    as209_median[.!outer_circle_mask] .= NaN

    # Save results
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cube.fits"), derotated...)
end

function rdi_subtract_epoch_4()
    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2005-07-27_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_2_aligned_frames.fits"), :)
    t222007 = load(joinpath(sequence_obslog.paths.sequences_folder, "t222007_4_aligned_frames.fits"), :)

    t222007_derotated = AstroImage[]
    for frame in t222007
        frame ./= frame["ITIME"]  # Normalize by integration time
        frame = rotate_image_center(frame, -frame["ROTPPOSN"])
        push!(t222007_derotated, copy(frame))
    end
    cube = framelist_to_cube(t222007)
    t222007_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

    search_radius = 5
    inner_mask_radius = 60
    outer_mask_radius = 230
    inner_circle_mask = make_circle_mask(size(as209[1]), inner_mask_radius)
    outer_circle_mask = make_circle_mask(size(as209[1]), outer_mask_radius)

    # Perform PSF subtraction with sub-pixel shifting, scaling, and offset
    derotated = AstroImage[]
    initial_guess = [0.0, 0.0, 1.0, 0.0]
    as209_derotated = AstroImage[]
    for frame in as209
        frame ./= frame["ITIME"]
        frame = rotate_image_center(frame, -frame["ROTPPOSN"])
        as209_derotated = push!(as209_derotated, copy(frame))

        @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."
        residual, optimal_params = optimal_subtract_target(frame, t222007_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

        residual = rotate_image_center(residual, residual["ROTPPOSN"])

        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle+residual["PARANG"])
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
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_median.fits"), as209_median)
    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_cube.fits"), derotated...)

    save(joinpath(sequence_obslog.paths.sequences_folder, "as209_derotated.fits"), as209_derotated...)
    save(joinpath(sequence_obslog.paths.sequences_folder, "t222007_derotated.fits"), t222007_derotated...)
end

@autolog begin

    rdi_subtract_epoch_1()
    rdi_subtract_epoch_2()
    rdi_subtract_epoch_3()
    rdi_subtract_epoch_4()

end