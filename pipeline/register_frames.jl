using Printf

using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations

using AIR
import AIR.crop 

function register_frames(frames, template_psf; sizes=(530, 515, 500))

    coarse_size, intermediate_size, final_size = sizes

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

    cropped_frames = AstroImage[]

    for frame in frames
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

    return fine_aligned_frames

end

function register_epoch_1()

    sequence_obslog_folder = "pipeline/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")
    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    unsaturated_sequences = ["as209_1", "hbc650_1", "hbc630_1"]
    for key in keys(sequence_obslog.sequences)

        if key in unsaturated_sequences
            @info "Skipping saturated sequence $key"
            continue
        end
        @info "Processing sequence" key

        frames = load_frames(sequence_obslog, key) 
        target = chop(key, tail=2)
        template_psf = load(joinpath(sequences_folder, "$(target)_1_template_psf_cored.fits"))

        sizes = (530, 515, 500)
        fine_aligned_frames = register_frames(frames, template_psf; sizes=sizes)
        save(joinpath(sequence_obslog.sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)

    end

end

@autolog begin
    register_epoch_1()
end