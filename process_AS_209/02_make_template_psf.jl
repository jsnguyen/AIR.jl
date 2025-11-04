using TOML
using AstroImages
using ImageTransformations
using Statistics
using CoordinateTransformations

using AIR
import AIR.crop 

function make_template_psf_and_center_frames(frames, sequences_folder; coarse_size=300, fine_size=250, kwargs...)

    cropped_frames, centered_frames = center_frames(frames, coarse_size, fine_size, first(frames)["ROTMODE"]; kwargs...)

    template_stack = framelist_to_cube(centered_frames)
    template_psf = median(template_stack, dims=3) |> x -> dropdims(x, dims=3)
    circle_mask = make_circle_mask(size(template_psf), 10)
    template_psf_cored = copy(template_psf)
    template_psf_cored[circle_mask] .= 0.0

    return cropped_frames, centered_frames, template_psf, template_psf_cored

end

function center_frames(frames, coarse_size, fine_size, rotator_mode; fixed_sigma=2.0, quantile_threshold=0.9999)

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
        cropped, _, _ = AIR.crop(frame, (coarse_size, coarse_size), center=coarse_center)
        push!(cropped_frames, cropped)
    end

    # center the frames
    centered_frames = Vector{AstroImage}(undef, length(cropped_frames))
    @info "Centering frames using Gaussian fitting..."

    Threads.@threads for i in eachindex(cropped_frames)
        frame = cropped_frames[i]

        # put frame into a "PSF-aligned frame" frame
        # vertical angle mode is already aligned
        if rotator_mode == "position angle"
            frame = rotate_image_center(frame, (frame["PARANG"]-frame["ROTPOSN"]))
            remove_nan!(frame)
        end

        max_idx = argquantile(frame.data, quantile_threshold)
        max_value = frame.data[max_idx...]
        initial_cy, initial_cx = Float64.(max_idx)
        
        inv_circle_mask = .! make_circle_mask(size(frame), 200; center=max_idx)
        background_level = median(frame.data[inv_circle_mask])
        
        initial_guess = [max_value, initial_cx, initial_cy, background_level]

        cropped_frame, final_cx, final_cy, _, _ = fit_and_crop(frame, (fine_size, fine_size), initial_guess; fixed_sigma=fixed_sigma)

        centered_frames[i] = cropped_frame
        @info "Frame $i" gaussian_center=(final_cy, final_cx)

        if final_cy < 0 || final_cx < 0 || final_cy > size(frame, 1) || final_cx > size(frame, 2)
            @warn "Frame $i center out of bounds, check fitting for frame $(frame["RED-FN"])"
        end
    end

    @info "Successfully centered $(length(centered_frames)) frames"

    return cropped_frames, centered_frames

end

@stage function make_template_psf(unsaturated_keys, sequences; coarse_size=400, fine_size=370, fixed_sigma=4.0, quantile_threshold=0.9999, use_cored=false)

    paths = context["paths"]

    template_psfs = Dict{String, AstroImage}()
    
    for key in unsaturated_keys
        @info "Making template PSF for sequence: $(length(sequences[key]))"
        cropped_frames, centered_frames, template_psf, template_psf_cored = make_template_psf_and_center_frames(sequences[key], paths.sequences_folder; coarse_size=coarse_size, fine_size=fine_size, fixed_sigma=fixed_sigma, quantile_threshold=quantile_threshold)

        remove_nan!(template_psf)
        max_idx = argquantile(template_psf, quantile_threshold)
        max_value = template_psf[max_idx...]
        initial_cy, initial_cx = Float64.(max_idx)

        inv_circle_mask = .! make_circle_mask(size(template_psf), Int(0.5*fine_size); center=max_idx)
        background_level = median(template_psf[inv_circle_mask])
        
        initial_guess = [max_value, initial_cx, initial_cy, fixed_sigma, fixed_sigma, background_level]

        res = fit_2d_gaussian(template_psf, initial_guess)
        amp, cx, cy, σx, σy, offset = res
        @info "Fitted template PSF parameters:" AMP=amp CX=cx CY=cy SIGMA_X=σx SIGMA_Y=σy OFFSET=offset

        template_psf = AstroImage(template_psf)

        template_psf["AMP"] = res[1]
        template_psf["CX"] = res[2]
        template_psf["CY"] = res[3]
        template_psf["SIGMA_X"] = res[4]
        template_psf["SIGMA_Y"] = res[5]
        template_psf["OFFSET"] = res[6]

        save(joinpath(paths.sequences_folder, "$(key)_cropped_sequence.fits"), cropped_frames...)
        save(joinpath(paths.sequences_folder, "$(key)_centered_sequence.fits"), centered_frames...)
        save(joinpath(paths.sequences_folder, "$(key)_template_psf.fits"), template_psf)
        save(joinpath(paths.sequences_folder, "$(key)_template_psf_cored.fits"), template_psf_cored)

        if use_cored
            template_psfs[key] = template_psf_cored
        else
            template_psfs[key] = template_psf
        end

    end

    return template_psfs

end