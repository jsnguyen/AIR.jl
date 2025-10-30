using Base.Threads
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Optim
using ADI
using Plots

using AIR
import AIR.crop 

function rdi_subtract(as209, reference, rotator_mode; search_radius=4.0, inner_mask_radius=60, outer_mask_radius=200)
    final_size = (2*outer_mask_radius, 2*outer_mask_radius)

    # derotate reference to PSF-aligned frame
    reference_psfup_cube = AstroImage[]
    for frame in reference
        frame,_,_ = crop(frame, final_size)
        frame ./= frame["ITIME"]  # Normalize by integration time
        if rotator_mode == "position angle"
            # rotate the reference to be in the pupil-up orientation
            # should have a diffraction spike aligned with the vertical
            frame = rotate_image_center(frame, (frame["PARANG"]-frame["ROTPOSN"]))
        end
        push!(reference_psfup_cube, copy(frame))
    end

    cube = framelist_to_cube(reference_psfup_cube)
    reference_psfup_median = median(cube, dims=3) |> x -> dropdims(x, dims=3)

    annulus_mask = make_annulus_mask(final_size, inner_mask_radius, outer_mask_radius)

    # Perform PSF subtraction with sub-pixel shifting, scaling, and offset
    as209_northup_cube = Vector{AstroImage}(undef, length(as209))
    as209_frameup_cube = Vector{AstroImage}(undef, length(as209))
    reference_frameup_cube = Vector{AstroImage}(undef, length(as209))
    initial_guess = [0.0, 0.0, 1.0, 0.0]
    Threads.@threads for i in eachindex(as209)
        frame = as209[i]
        frame,_,_ = crop(frame, final_size)

        frame ./= frame["ITIME"]

        # rotate the reference to be in the as209 frame orientation
        # should be aligned with the as209 frame
        if rotator_mode == "position angle"
            rerotated_reference_psfup_median = rotate_image_center(reference_psfup_median, -(frame["PARANG"]-frame["ROTPOSN"]))
        else
            # dont need to rotate in vertical angle mode
            rerotated_reference_psfup_median = reference_psfup_median
        end

        as209_frameup_cube[i] = deepcopy(frame)
        reference_frameup_cube[i] = deepcopy(rerotated_reference_psfup_median)

        @info "Performing PSF subtraction with sub-pixel shift, scale, and offset optimization..."
        residual, optimal_params = optimal_subtract_target(frame, rerotated_reference_psfup_median, initial_guess, search_radius; inner_mask_radius=inner_mask_radius, outer_mask_radius=outer_mask_radius)

        angle, _ = calculate_north_angle(residual.header)
        derot = rotate_image_center(residual, -angle)

        derot[.!annulus_mask] .= NaN
        as209_northup_cube[i] = derot

        initial_guess = optimal_params

        @info "Optimal Params" optimal_params=optimal_params
    end

    as209_northup_median = framelist_to_cube(as209_northup_cube) |> x -> median(x, dims=3) |> x -> dropdims(x, dims=3)
    as209_northup_median[.!annulus_mask] .= NaN

    return as209_northup_median, as209_northup_cube, as209_frameup_cube, reference_psfup_median, reference_psfup_cube, reference_frameup_cube

end

function rdi_subtract_position_angle_epochs()

    for date in ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]
        sequence_obslog_path = joinpath("pipeline/obslogs/$(date)_sequences.toml")

        @info "Loading sequence_obslog from" sequence_obslog_path
        sequence_obslog = Obslog(sequence_obslog_path)

        if date == "2002-06-16"
            rotator_mode = "vertical angle"
            as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_4_aligned_frames.fits"), :)
            reference = load(joinpath(sequence_obslog.paths.sequences_folder, "hbc650_3_aligned_frames.fits"), :)
        elseif date == "2002-08-02"
            rotator_mode = "position angle"
            as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_2_aligned_frames.fits"), :)
            reference = deepcopy(as209) # enough movement in frame to use itself as reference
        elseif date == "2002-08-21"
            rotator_mode = "position angle"
            as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_4_aligned_frames.fits"), :)
            reference = load(joinpath(sequence_obslog.paths.sequences_folder, "tyc2307_3_aligned_frames.fits"), :)
        elseif date == "2005-07-27"
            rotator_mode = "position angle"
            as209 = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_2_aligned_frames.fits"), :)
            reference = load(joinpath(sequence_obslog.paths.sequences_folder, "t222007_4_aligned_frames.fits"), :)
        end

        as209_northup_median, as209_northup_cube, as209_frameup_cube, reference_psfup_median, reference_psfup_cube, reference_frameup_cube  = rdi_subtract(as209, reference, rotator_mode)

        # Save results
        # northup -> north up on frame
        # psfup -> PSF diffraction spike aligned with the vertical
        # frameup -> aligned with the as209 frame, reference rotated to match

        # final data products
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_median.fits"), as209_northup_median)
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_cube.fits"), as209_northup_cube...)

        # median PSF-aligned frame we are using for subtracting
        # psfup means that the parallactic angle and the user set position removed so that we are aligned with the pupil
        # basically turning it back into vertical angle mode
        save(joinpath(sequence_obslog.paths.sequences_folder, "reference_psfup_median.fits"), reference_psfup_median)
        save(joinpath(sequence_obslog.paths.sequences_folder, "reference_psfup_cube.fits"), reference_psfup_cube...)

        # these are the aligned frames we are actually minimizing residuals against
        # frameup means the reference is rotated to match the target frame
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_frameup_cube.fits"), as209_frameup_cube...)
        save(joinpath(sequence_obslog.paths.sequences_folder, "reference_frameup_cube.fits"), reference_frameup_cube...)
    end

end

@autolog begin

    rdi_subtract_position_angle_epochs()

end