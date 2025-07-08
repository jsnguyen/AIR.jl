using AstroImages
using ImageTransformations
using Statistics
using CoordinateTransformations

using AIR

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences = load_sequences(sequence_obslog)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")

    unsaturated_sequences = ["as209_1", "hbc650_1", "hbc630_1"]

    coarse_size = 300
    fine_size = 256

    for key in keys(sequences)
        if !(key in unsaturated_sequences)
            @info "Skipping saturated sequence $key"
            continue
        end


        # subtract off mean background

        mean_background = 0.0
        for frame in sequences[key]
            mean_background += measure_background(frame; mask_radius=200)
        end
        mean_background /= length(sequences[key])

        for frame in sequences[key]
            frame .-= mean_background
            frame["BGSUB"] = mean_background
        end

        # crop to speed things up a bit
        cropped_frames = AstroImage[]
        for frame in sequences[key]
            _, coarse_center = findmax(frame.data)
            cropped, _, _ = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
            push!(cropped_frames, cropped)
        end

        # Center each frame individually using Gaussian fitting
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
        
        # Create template PSF from centered frames
        if length(centered_frames) > 0
            template_stack = framelist_to_cube(centered_frames)
            template_psf = median(template_stack, dims=3) |> x->dropdims(x, dims=3)
            circle_mask = make_circle_mask(size(template_psf), 10)
            template_psf_cored = copy(template_psf)  # Create a copy to modify
            template_psf_cored[circle_mask] .= 0.0

            
            # Save results
            save(joinpath(sequences_folder, "$(key)_cropped_sequence.fits"), cropped_frames...)
            save(joinpath(sequences_folder, "$(key)_centered_sequence.fits"), centered_frames...)
            save(joinpath(sequences_folder, "$(key)_template_psf.fits"), template_psf)
            save(joinpath(sequences_folder, "$(key)_template_psf_cored.fits"), template_psf_cored)
        end

    end

end