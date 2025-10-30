using Printf

using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations

using AIR
import AIR.crop 

function register_frames(frames, template_psf, rotator_mode; sizes=(530, 515, 500))

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
        coarse_center = argquantile(frame.data, 0.99999)
        cropped,_,_ = crop(frame, (coarse_size, coarse_size); center=Tuple(coarse_center))
        push!(cropped_frames, cropped)
    end

    @info "Mean background: $(@sprintf("%.3f", mean_background))"

    n_sigma_ccr = 30.0

    aligned_frames = AstroImage[]

    @showprogress desc="Align" for (i,frame) in enumerate(cropped_frames)

        # rotate the template PSF to match the frame's rotation
        if rotator_mode == "position angle"
            template_psf = rotate_image_center(template_psf, -(frame["PARANG"]-frame["ROTPOSN"]))
            remove_nan!(template_psf)
        end

        aligned = cross_correlate_align(frame, template_psf, n_sigma_ccr)

        aligned,_ ,_ =  crop(aligned, (intermediate_size, intermediate_size))
        remove_nan!(aligned)

        push!(aligned_frames, aligned)
    end

    fine_aligned_frames = AstroImage[]
    @showprogress desc="Fine align" for frame in aligned_frames

        aligned = cross_correlate_align(frame, aligned_frames[1], n_sigma_ccr)

        aligned,_ ,_ =  crop(aligned, (final_size, final_size))
        remove_nan!(aligned)

        push!(fine_aligned_frames, aligned)
    end

    return fine_aligned_frames

end

function register_epochs()

    dates = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

    picked_sequences  = [["as209_2", "as209_3", "as209_4", "hbc630_2", "hbc650_2", "hbc650_3"],
                         ["as209_2"],
                         ["as209_3", "as209_4", "tyc2307_2", "tyc2307_3"],
                         ["as209_2", "t222007_4"]]

    registering_sizes = ((530, 500, 400),
                         (450, 430, 400),
                         (420, 410, 400),
                         (420, 410, 400))

    sequence_obslog_folder = "pipeline/obslogs"
    for (sizes,picked_sequence,date) in zip(registering_sizes, picked_sequences, dates)
        sequence_obslog_path = joinpath(sequence_obslog_folder, "$(date)_sequences.toml")
        @info "Loading sequence_obslog from" sequence_obslog_path
        sequence_obslog = Obslog(sequence_obslog_path)

        for key in keys(sequence_obslog.sequences)

            if !(key in picked_sequence)
                @info "Skipping sequence $key"
                continue
            end
            @info "Processing sequence" key

            if date=="2002-08-02"
                template_psf = load("data/2002-06-16/sequences/as209_template_psf_cored.fits")
            else
                target = chop(key, tail=2)
                template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "$(target)_template_psf_cored.fits"))
            end

            frames = load_frames(sequence_obslog, key) 
            fine_aligned_frames = register_frames(frames, template_psf, frames[1]["ROTMODE"]; sizes=sizes)
            save(joinpath(sequence_obslog.paths.sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)

        end

    end

end

@autolog begin
    register_epochs()
end