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

function register_frames(template_psf, frames; sizes=(200, 190, 180))

    rotator_mode = first(frames)["ROTMODE"]

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
        median_frame = mapwindow(x -> median(vec(x)), frame.data, (3,3))
        coarse_center = argquantile(median_frame, 0.99999)
        cropped,_,_ = crop(frame, (coarse_size, coarse_size); center=Tuple(coarse_center))
        @info "Centering frame $(frame["RED-FN"]) at $coarse_center"
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

@stage function register_sequences(template_psf_key, target_keys, sequences, template_psfs; injected_psf=nothing, kwargs...)

    paths = context["paths"]

    if injected_psf !== nothing
        @info "Using injected PSF for registration"
        template_psf = injected_psf
    else
        @info "Using template PSF key: $template_psf_key"
        template_psf = template_psfs[template_psf_key]
    end

    registered_frames = Dict{String, Any}()
    for key in target_keys
        frames = sequences[key]
        rf = register_frames(template_psf, frames; kwargs...)
        save(joinpath(paths.sequences_folder, "$(key)_aligned_frames.fits"), rf...)
        registered_frames[key] = rf
    end
    
    return registered_frames

end