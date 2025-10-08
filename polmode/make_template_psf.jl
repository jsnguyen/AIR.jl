using TOML
using AstroImages
using ImageTransformations
using Statistics
using CoordinateTransformations

using AIR

function make_template_psf_from_sequence(sequence, key, sequences_folder; coarse_size=300, fine_size=256, unsaturated_sequences=String[], kwargs...)

    cropped_frames, centered_frames = make_template_psf(sequence, coarse_size, fine_size, sequence[1]["ROTMODE"]; kwargs...)

    # make final template
    if length(centered_frames) > 0
        template_stack = framelist_to_cube(centered_frames)
        template_psf = median(template_stack, dims=3) |> x -> dropdims(x, dims=3)
        circle_mask = make_circle_mask(size(template_psf), 10)
        template_psf_cored = copy(template_psf)
        template_psf_cored[circle_mask] .= 0.0

        key = chop(key, tail=2)  # remove the _1 or _2 suffix
        save(joinpath(sequences_folder, "$(key)_cropped_sequence.fits"), cropped_frames...)
        save(joinpath(sequences_folder, "$(key)_centered_sequence.fits"), centered_frames...)
        save(joinpath(sequences_folder, "$(key)_template_psf.fits"), template_psf)
        save(joinpath(sequences_folder, "$(key)_template_psf_cored.fits"), template_psf_cored)
    end

end

function make_template_psf(frames, coarse_size, fine_size, rotator_mode; fixed_sigma=2.0, quantile_threshold=0.9999)

    # subtract off mean background
    mean_background = 0.0
    for frame in frames
        mean_background += measure_background(frame; mask_radius=200)
    end
    mean_background /= length(frames)

    for frame in frames
        frame .-= mean_background
        frame["BGSUB"] = mean_background
    end

    # crop to speed things up a bit
    cropped_frames = AstroImage[]
    for frame in frames
        @info "Processing" frame["RED-FN"]
        coarse_center = argquantile(frame.data, quantile_threshold)
        cropped, _, _ = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
        push!(cropped_frames, cropped)
    end

    # center the frames
    centered_frames = AstroImage[]
    @info "Centering frames using Gaussian fitting..."

    for (i, frame) in enumerate(cropped_frames)

        # put frame into a "PSF-aligned frame" frame
        # vertical angle mode is already aligned
        if rotator_mode == "position angle"
            frame = rotate_image_center(frame, (frame["PARANG"]-frame["ROTPOSN"]))
            remove_nan!(frame)
        end

        max_idx = argquantile(frame.data, quantile_threshold)
        initial_cy, initial_cx = Float64.(Tuple(max_idx))
        initial_guess = [5000.0, initial_cx, initial_cy, 10.0]

        cropped_frame, final_cx, final_cy, _, _ = fit_and_crop(frame, (fine_size, fine_size), initial_guess; fixed_sigma=fixed_sigma)

        push!(centered_frames, cropped_frame)
        @info "Frame $i" gaussian_center=(final_cy, final_cx)
    end

    @info "Successfully centered $(length(centered_frames)) frames"

    return cropped_frames, centered_frames

end

function make_template_psf_epochs()
    
    date = "2025-10-07"
    coarse_size = 220
    fine_size = 120

    sequence_obslog_path = "live_ingestion/obslogs/$(date)_split_sequences.toml"
    @info "Loading sequence_obslog from" sequence_obslog_path

    sequence_obslog = Obslog(sequence_obslog_path)

    key = "HD215806_dn_1"
    make_template_psf_from_sequence(sequence_obslog.sequences[key], key, sequence_obslog.paths.sequences_folder; coarse_size=coarse_size, fine_size=fine_size)

end

@autolog begin
    make_template_psf_epochs()
end