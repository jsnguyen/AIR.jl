using TOML
using AstroImages
using ImageTransformations
using Statistics
using CoordinateTransformations

using AIR

function make_template_psf_from_sequences(sequence_obslog; coarse_size=300, fine_size=256, unsaturated_sequences=String[])

    sequences = sequence_obslog.sequences

    for key in keys(sequences)
        if !(key in unsaturated_sequences)
            @info "Skipping saturated sequence $key"
            continue
        end

        cropped_frames, centered_frames = make_template_psf(sequences[key], coarse_size, fine_size)

        # make final template
        if length(centered_frames) > 0
            template_stack = framelist_to_cube(centered_frames)
            template_psf = median(template_stack, dims=3) |> x -> dropdims(x, dims=3)
            circle_mask = make_circle_mask(size(template_psf), 10)
            template_psf_cored = copy(template_psf)
            template_psf_cored[circle_mask] .= 0.0

            save(joinpath(sequence_obslog.sequences_folder, "$(key)_cropped_sequence.fits"), cropped_frames...)
            save(joinpath(sequence_obslog.sequences_folder, "$(key)_centered_sequence.fits"), centered_frames...)
            save(joinpath(sequence_obslog.sequences_folder, "$(key)_template_psf.fits"), template_psf)
            save(joinpath(sequence_obslog.sequences_folder, "$(key)_template_psf_cored.fits"), template_psf_cored)
        end

    end

end

function make_template_psf(frames, coarse_size, fine_size)

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
        _, coarse_center = findmax(frame.data)
        cropped, _, _ = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
        push!(cropped_frames, cropped)
    end

    # center the frames
    centered_frames = AstroImage[]
    @info "Centering frames using Gaussian fitting..."

    for (i, frame) in enumerate(cropped_frames)

        _, max_idx = findmax(frame)
        initial_cy, initial_cx = Float64.(Tuple(max_idx))
        initial_guess = [5000.0, initial_cx, initial_cy, 10.0]

        cropped_frame, final_cx, final_cy, _, _ = fit_and_crop(frame.data, (fine_size, fine_size), initial_guess; fixed_sigma=2.0)

        push!(centered_frames, AstroImage(cropped_frame, frame.header))
        @info "Frame $i" gaussian_center=(final_cy, final_cx)
    end

    @info "Successfully centered $(length(centered_frames)) frames"

    return cropped_frames, centered_frames

end

function make_template_psf_epoch_1()

    obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(obslog_folder, "2002-06-16_sequences.toml")
    rejects_obslog_path = joinpath(obslog_folder, "2002-06-16_rejects.toml")
    rejects = load_rejects(rejects_obslog_path)

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path, rejects=rejects)

    unsaturated_sequences = ["as209_1", "hbc650_1", "hbc630_1"]

    coarse_size = 300
    fine_size = 256

    make_template_psf_from_sequences(sequence_obslog; coarse_size=coarse_size, fine_size=fine_size, unsaturated_sequences=unsaturated_sequences)

end

@autolog begin
    make_template_psf_epoch_1()
end