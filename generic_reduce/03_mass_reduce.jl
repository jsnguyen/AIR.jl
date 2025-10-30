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

    @inbounds for index in bad_indices
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
function reduce_frame(frame, master_flats, master_darks, master_skies, masks;median_size = 6, gain=1.0)

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

    if length(master_skies) > 0
        ind, matched_sky = find_closest_sky(frame, master_skies)
        if matched_sky !== nothing
            @info "Finding sky frame for $(reduced["FILENAME"])"
            @info "Subtracting sky frame from $(reduced["FILENAME"])"
            reduced .-= matched_sky
            reduced["SKYSUB"] = true
        else
            reduced["SKYSUB"] = false
        end
    else
        reduced["SKYSUB"] = false
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

@stage function mass_reduce_frames(master_darks, master_flats, master_skies, master_masks; guest_obslog=nothing, use_guest_flats=false, use_guest_darks=false, use_guest_skies=false, use_guest_masks=false)

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

    reduced_frames = Vector{AstroImage}(undef, length(sci))
    Threads.@threads for i in eachindex(sci)

        frame = sci[i]

        # "guest" master frames will supersede the main obslog's master frames if specified
        master_flats = use_guest_flats ? vcat(guest_obslog.master_flats, master_flats) : master_flats
        master_darks = use_guest_darks ? vcat(guest_obslog.master_darks, master_darks) : master_darks
        master_skies = use_guest_skies ? vcat(guest_obslog.master_skies, master_skies) : master_skies
        master_masks = use_guest_masks ? vcat(guest_obslog.master_masks, master_masks) : master_masks

        reduced = reduce_frame(frame, master_flats, master_darks, master_skies, master_masks; gain=gain)

        readnoise = get_NIRC2_readnoise(reduced["SAMPMODE"])
        reduced["RDNOISE"] = readnoise

        reduced_frames[i] = reduced

    end

    @info "Reduced $(length(reduced_frames)) science frames"

    make_and_clear(paths.reduced_folder, "reduced_*.fits")
    
    @info "Saving reduced frames..."
    for rf in reduced_frames
        reduced_filepath = joinpath(paths.reduced_folder, rf["RED-FN"])
        @context_save key="reduced_frames" save(reduced_filepath, rf)
    end

    return reduced_frames

end
