using Printf

using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations

using AIR
import AIR.crop 

@autolog begin

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")
    sequences = load_sequences(sequence_obslog)

    unsaturated_sequences = ["as209_1", "hbc650_1", "hbc630_1"]
    for key in keys(sequences)

        if key in unsaturated_sequences
            @info "Skipping saturated sequence $key"
            continue
        end
        @info "Processing sequence" key

        frames = load_frames(sequence_obslog, key) 
        target = chop(key, tail=2)
        template_psf = load(joinpath(sequences_folder, "$(target)_1_template_psf_cored.fits"))

        coarse_size = 530 # size of the coarse crop
        intermediate_size = 515 # size of the intermediate crop
        final_size = 500

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

        cropped_frames = AstroImage[]

        for frame in sequences[key]
            _, coarse_center = findmax(frame.data)
            cropped,_,_ = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
            push!(cropped_frames, cropped)
        end

        @info "Mean background: $(@sprintf("%.3f", mean_background))"

        n_sigma_ccr = 30.0

        aligned_frames = AstroImage[]

        @showprogress desc="Align" for (i,frame) in enumerate(cropped_frames)

            aligned = cross_correlate_align(frame.data, template_psf, n_sigma_ccr)

            aligned,_ ,_ =  crop(aligned, (intermediate_size, intermediate_size))
            @. aligned[~isfinite(aligned)] = 0
            aligned = AstroImage(aligned, frame.header)

            push!(aligned_frames, aligned)
        end

        fine_aligned_frames = AstroImage[]
        @showprogress desc="Fine align" for frame in aligned_frames

            aligned = cross_correlate_align(frame.data, aligned_frames[1], n_sigma_ccr)

            aligned,_ ,_ =  crop(aligned, (final_size, final_size))
            @. aligned[~isfinite(aligned)] = 0
            aligned = AstroImage(aligned, frame.header)

            push!(fine_aligned_frames, aligned)
        end

        save(joinpath(sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)

    end

end