using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics
using ProgressMeter

using AIR

```
median replace using local pixels
if no good pixels, replace with fail_val
```
function local_median_replace_bad_pixels!(data, mask, median_size; fail_val=0.0)
    bad_indices = findall(mask)
    half_size = median_size ÷ 2
    height, width = size(data)

    @inbounds for idx in bad_indices
        i, j = Tuple(idx)

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
function reduce_frame(frame, master_flats, master_darks, masks; median_size = 6, gain=1.0)

    reduced = copy(frame)

    matched_flat = find_closest_flat(reduced, master_flats)
    matched_dark = find_closest_dark(reduced, master_darks)

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

@autolog begin

    @info "Using $(Threads.nthreads()) threads for reduction"

    allowed = ["2025-10-07"]

    obslog_folder = "polmode/obslogs"
    for obslog_filename in Glob.glob("*_obslog.toml", obslog_folder)

        if !any(date -> occursin(date, obslog_filename), allowed)
            @warn "Skipping file: $obslog_filename (not in allowed dates)"
            continue
        end

        @info "Loading obslog from" obslog_filename
        obslog = Obslog(obslog_filename)

        p = Progress(length(obslog.sci))
        update!(p, 0)
        counter = Threads.Atomic{Int}(0)
        lock = Threads.SpinLock()

        gain = get_NIRC2_gain(obslog.date) # gain is the same for all dates

        # no LACosmic, probably more trouble than it's worth

        reduced_frames = Vector{AstroImage}(undef, length(obslog.sci))
        Threads.@threads for i in eachindex(obslog.sci)

            frame = obslog.sci[i]

            reduced = reduce_frame(frame, obslog.master_flats, obslog.master_darks, obslog.masks; gain=gain)

            readnoise = get_NIRC2_readnoise(reduced["SAMPMODE"])
            reduced["RDNOISE"] = readnoise

            reduced_frames[i] = reduced

            Threads.atomic_add!(counter, 1)
            Threads.lock(lock)
            update!(p, counter[])
            Threads.unlock(lock)
        end

        @info "Reduced $(length(reduced_frames)) science frames"

        make_and_clear(obslog.paths.reduced_folder, "reduced_*.fits")

        reduced_filepaths = String[]
        for rf in reduced_frames
            reduced_filepath = joinpath(obslog.paths.reduced_folder, rf["RED-FN"])
            push!(reduced_filepaths, rf["RED-FN"])
            @info "Saving reduced frame to $(reduced_filepath)"
            save(reduced_filepath, rf)
        end

        reduced_obslog = OrderedDict{String, Any}("data_folder" => obslog.paths.data_folder,
                                                  "date" => obslog.date,
                                                  "reduced" => OrderedDict{String, Any}(
                                                  "reduced_sci" => reduced_filepaths))

        write_toml(obslog.paths.reduced_file, reduced_obslog)

    end

end