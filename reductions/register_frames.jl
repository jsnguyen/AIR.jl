using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LsqFit

using AIR
import AIR.crop 

function align_to_template(frame::AstroImage, template::AstroImage; σ::Real=5.0, fillval=NaN)
    # 1) cross‐correlate
    cc = imfilter(frame.data, centered(template.data))

    # 2) find integer peak
    _, ci = findmax(cc)
    r, c = Tuple(ci)

    # 3) sub‐pixel Gaussian fit: p = [A, x0, y0, σ]
    p0  = [cc[ci], Float64(c), Float64(r), σ]
    fit = gaussian2d_fixedwidth_fit(cc, p0, σ)
    x0, y0 = fit[2], fit[3]

    # 4) compute 1‐based center of cc

    # no idea why there's a -2 here but it works better
    # something to do with indexing and convolution sizes?
    cy = (size(cc,1)-2) / 2
    cx = (size(cc,2)-2) / 2

    # 5) offset = (Δrow, Δcol)
    #offset = (y0 - cy, x0 - cx)
    offset = (y0 - cy, x0 - cx)
    @info "offset" offset=offset

    # 6) warp the raw data by this sub‐pixel shift
    warped = warp(frame.data,
                  Translation(offset...),
                  axes(frame.data),
                  fill=fillval)

    # 7) return with original header
    return AstroImage(warped, frame.header)
end

"""
    measure_background(frame::AstroImage; mask_radius=50, edge_buffer=20)

Measure the background level in a frame while masking out the PSF.
Returns the median background level from an annular region.
"""
function measure_background(frame::AstroImage; mask_radius=50)
    data = frame.data
    rows, cols = size(data)
    
    # Find the center of the PSF (brightest pixel)
    _, center_idx = findmax(data)
    cy, cx = Tuple(center_idx)
    
    # Create mask to exclude PSF and edges
    mask = trues(size(data))
    
    # Mask out the PSF (circular region around center)
    for i in 1:rows, j in 1:cols
        r = sqrt((i - cy)^2 + (j - cx)^2)
        if r < mask_radius
            mask[i, j] = false
        end
    end
    
    # Extract background pixels
    background_pixels = data[mask]
    
    if length(background_pixels) == 0
        @warn "No background pixels found, returning 0"
        return 0.0
    end
    
    # Return median background level
    return median(background_pixels)
end

autolog("$(@__FILE__).log") do

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
        
        # Arrays to store background measurements
        background_levels = Float64[]
        
        # Align the frames to the template PSF

        cropped_frames = AstroImage[]
        for frame in sequences[key]
            _, coarse_center = findmax(frame.data)
            cropped = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
            push!(cropped_frames, cropped)
            
            # Measure background level in the cropped frame
            bg_level = measure_background(cropped; mask_radius=200)
            push!(background_levels, bg_level)
        end

        mean_background = mean(background_levels)
        
        @info "Background levels for sequence $key: $(round.(background_levels, digits=3))"
        @info "Mean background: $(round(mean_background, digits=3)), Std: $(round(std(background_levels), digits=3))"

        n_sigma_ccr = 50.0

        aligned_frames = AstroImage[]
        for frame in cropped_frames
            aligned = align_to_template(frame, template_psf; σ=n_sigma_ccr, fillval=0.0)
            aligned = aligned |> x -> crop(x, (intermediate_size, intermediate_size))
            # Replace NaN or Inf values with zero
            aligned.data .= ifelse.(isfinite.(aligned.data), aligned.data, 0.0)
            push!(aligned_frames, aligned)
        end

        fine_aligned_frames = AstroImage[]
        for frame in aligned_frames
            aligned = align_to_template(frame, aligned_frames[1]; σ=n_sigma_ccr, fillval=0.0)
            aligned = aligned |> x -> crop(x, (final_size, final_size)) |> x -> x.-mean_background
            # Replace NaN or Inf values with zero
            aligned.data .= ifelse.(isfinite.(aligned.data), aligned.data, 0.0)
            push!(fine_aligned_frames, aligned)
        end

        save(joinpath(sequences_folder, "$(key)_aligned_frames.fits"), fine_aligned_frames...)
        
    end

end