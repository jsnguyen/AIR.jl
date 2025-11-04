using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics
using Optim

using AIR

#median replace using local pixels
#if no good pixels, replace with fail_val
function local_median_replace_bad_pixels!(data, mask, median_size; fail_val=0.0)
    bad_indices = findall(mask)
    half_size = median_size ÷ 2
    height, width = size(data)

    for index in bad_indices
        i, j = Tuple(index)

        i_start = max(1, i - half_size)
        i_end = min(height, i + half_size)
        j_start = max(1, j - half_size)
        j_end = min(width, j + half_size)

        window = @view data[i_start:i_end, j_start:j_end]
        window_mask = @view mask[i_start:i_end, j_start:j_end]

        # Fast median of good pixels
        good_pixels = window[.!window_mask]
        if length(good_pixels) > 0
            data[i, j] = median(good_pixels)
        else
            data[i, j] = fail_val
        end
    end

end

# median_size is the median pixel replacement size, even numbers work best
function reduce_frame(frame, master_flats, master_darks, master_skies, masks;median_size = 6, gain=1.0, skip_sky_sub=false)

    reduced = copy(frame)

    flat_ind, matched_flat = find_closest_flat(reduced, master_flats)
    dark_ind, matched_dark = find_closest_dark(reduced, master_darks)

    if matched_dark === nothing
        #@warn "No matching dark found for $(reduced["FILENAME"])" reduced["ITIME"] reduced["COADDS"]
        reduced["DARKSUB"] = false
        matched_dark = zeros(size(reduced))
    else
        reduced["DARKSUB"] = true
    end

    if matched_flat === nothing
        #@warn "No matching flat found for $(reduced["FILENAME"])" reduced["FILTER"]
        reduced["FLATDIV"] = false
        matched_flat = ones(size(reduced))
    else
        reduced["FLATDIV"] = true
    end

    reduced["SKYSUB"] = false
    if length(master_skies) > 0 && !skip_sky_sub
        ind, matched_sky = find_closest_sky(frame, master_skies)
        if matched_sky !== nothing
            @info "Finding sky frame for $(reduced["FILENAME"])"
            @info "Subtracting sky frame from $(reduced["FILENAME"])"
            reduced .-= matched_sky
            reduced["SKYSUB"] = true
        end
    end

    reduced = (reduced .- matched_dark) ./ matched_flat

    reduced ./= reduced["COADDS"] # divide by coadds
    reduced .*= gain

    # start with the bad pixel mask as our mask
    mask = copy(NIRC2_bad_pixel_mask)

    if size(mask) != size(reduced)
        mask, _, _ = crop(NIRC2_bad_pixel_mask, size(reduced))
    end

    if haskey(masks, size(reduced))
        mask = mask .| masks[size(reduced)] # also combine the extra mask if it exists
    end

    nan_mask = isnan.(reduced.data) .| isinf.(reduced.data) # create a mask for NaN and Inf values
    mask = mask .| nan_mask # combine the bad pixel mask with the NaN/Inf mask

    local_median_replace_bad_pixels!(reduced.data, mask, median_size)

    reduced_filename = "reduced_$(lpad(reduced["FRAMENO"], 4, '0')).fits"
    reduced["RED-FN"] = reduced_filename

    return reduced

end

@stage function mass_reduce_frames(master_darks, master_flats, master_skies, master_masks; guest_master_flats=nothing, guest_master_darks=nothing, guest_master_skies=nothing, kwargs...)

    date = context["paths"].date
    paths = context["paths"]
    sci_filepaths = register["01"][1]

    if isfile(context["paths"].rejects_file)
        rejects = load_rejects(context["paths"].rejects_file)
    else
        rejects = String[]
    end

    sci = load_frames(sci_filepaths; rejects=rejects)
    @info "Loaded $(length(sci)) science frames for reduction"
    @info "Rejects $(length(rejects)) frames based on rejects file"

    gain = get_NIRC2_gain(date) # gain is the same for all dates

    # no LACosmic, probably more trouble than it's worth

    # "guest" master frames will supersede the main obslog's master frames if specified
    master_flats = guest_master_flats !== nothing ? vcat(guest_master_flats, master_flats) : master_flats
    master_darks = guest_master_darks !== nothing ? vcat(guest_master_darks, master_darks) : master_darks
    master_skies = guest_master_skies !== nothing ? vcat(guest_master_skies, master_skies) : master_skies

    reduced_frames = Vector{AstroImage}(undef, length(sci))
    Threads.@threads for i in eachindex(sci)

        frame = sci[i]

        reduced = reduce_frame(frame, master_flats, master_darks, master_skies, master_masks; gain=gain, kwargs...)

        readnoise = get_NIRC2_readnoise(reduced["SAMPMODE"])
        reduced["RDNOISE"] = readnoise

        reduced_frames[i] = reduced

    end

    @info "Reduced $(length(reduced_frames)) science frames"

    make_and_clear(paths.reduced_folder, "reduced_*.fits")
    
    @info "Saving reduced frames..."

    reduced_filepaths = [joinpath(paths.reduced_folder, rf["RED-FN"]) for rf in reduced_frames]
    @context_store context reduced_filepaths
    Threads.@threads for i in eachindex(reduced_frames)
        fp = reduced_filepaths[i]
        rf = reduced_frames[i]
        save(fp, rf)
    end

    return reduced_frames

end
