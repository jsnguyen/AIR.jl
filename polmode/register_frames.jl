using Printf

using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Serialization

using AIR
import AIR.crop 

function register_frames(frames, template_psf, rotator_mode; sizes=(240, 230, 220, 210))

    coarse_size, intermediate_size, final_size = sizes

    # subtract off mean background
    mean_background = 0.0
    for frame in frames
        mean_background += measure_background(frame; mask_radius=100)
    end
    mean_background /= length(frames)

    for frame in frames
        frame .-= mean_background
        frame["BGSUB"] = mean_background
    end

    cropped_frames = AstroImage[]

    for frame in frames
        coarse_center = argquantile(frame.data, 0.99)
        cropped,_,_ = crop(frame, (coarse_size, coarse_size); center=Tuple(coarse_center))
        push!(cropped_frames, cropped)
    end

    @info "Mean background: $(@sprintf("%.3f", mean_background))"

    n_sigma_ccr = 10.0

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

    hyper_fine_aligned_frames = AstroImage[]
    @showprogress desc="Hyper Fine align" for frame in aligned_frames

        circle_mask = make_circle_mask(size(fine_aligned_frames[1]), 10)
        cored_frame = copy(fine_aligned_frames[1])
        cored_frame[circle_mask] .= 0.0

        aligned = cross_correlate_align(frame, cored_frame, n_sigma_ccr)

        aligned,_ ,_ =  crop(aligned, (final_size, final_size))
        remove_nan!(aligned)

        push!(hyper_fine_aligned_frames, aligned)
    end

    angles = Float64[]
    for frame in frames
        ang , _ = calculate_north_angle(frame.header)
        push!(angles, ang)
    end

    return hyper_fine_aligned_frames, angles

end

function register_epochs()

    date = "2025-10-07"

    sequence_obslog_folder = "live_ingestion/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "$(date)_split_sequences.toml")
    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    key = "BD45598_up_1"

    @info "Processing sequence" key

    #target = chop(key, tail=2)
    #template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "$(target)_template_psf.fits"))
    template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "HD215806_up_template_psf_cored.fits"))

    frames = load_frames(sequence_obslog, key) 
    fine_aligned_frames, angles = register_frames(frames, template_psf, frames[1]["ROTMODE"])
    @info "ANGLES" angles

    open(joinpath(sequence_obslog.paths.sequences_folder, "angles.jls"), "w") do io
        serialize(io, angles)
    end

    save(joinpath(sequence_obslog.paths.sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)

    key = "BD45598_dn_1"

    @info "Processing sequence" key

    #target = chop(key, tail=2)
    #template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "$(target)_template_psf.fits"))
    template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "HD215806_up_template_psf_cored.fits"))

    frames = load_frames(sequence_obslog, key) 
    fine_aligned_frames, angles = register_frames(frames, template_psf, frames[1]["ROTMODE"])
    save(joinpath(sequence_obslog.paths.sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)


end

@autolog begin
    register_epochs()
end