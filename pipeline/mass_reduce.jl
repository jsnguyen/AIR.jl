using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics

using AIR

@autolog begin

    median_size = 5 # median pixel replacement size

    @info "Using $(Threads.nthreads()) threads for reduction"

    obslog_folder = "pipeline/obslogs"
    for obslog_filename in Glob.glob("*_obslog.toml", obslog_folder)
        @info "Loading obslog from" obslog_filename

        obslog = Obslog(obslog_filename)
        master_darks = obslog.master_darks
        master_flats = obslog.master_flats
        masks = obslog.masks
        sci = obslog.sci

        # add cosmic ray image cleaning here at some point?
        # first try didn't work properly...

        reduced_frames = Vector{AstroImage}(undef, length(sci))
        Threads.@threads for i in eachindex(sci)
            sf = sci[i]

            matched_flat = find_closest_flat(sf, master_flats)
            matched_dark = find_closest_dark(sf, master_darks)
            
            if matched_dark === nothing
                @warn "No matching dark found for $(sf["FILENAME"])" sf["ITIME"] sf["COADDS"]
                sf["DARKSUB"] = false
                matched_dark = zeros(size(sf))
            else
                sf["DARKSUB"] = true
            end

            if matched_flat === nothing
                @warn "No matching flat found for $(sf["FILENAME"])" sf["FILTER"]
                sf["FLATDIV"] = false
                matched_flat = ones(size(sf))
            else
                sf["FLATDIV"] = true
            end

            # start with the bad pixel mask as our mask
            mask = copy(NIRC2_bad_pixel_mask)

            if size(mask) != size(sf)
                mask, _, _ = crop(NIRC2_bad_pixel_mask, size(sf))
            end

            if haskey(masks, size(sf))
                mask = mask .| masks[size(sf)] # also combine the extra mask if it exists
            end

            nan_mask = isnan.(sf.data) .| isinf.(sf.data) # create a mask for NaN and Inf values
            mask = mask .| nan_mask # combine the bad pixel mask with the NaN/Inf mask

            # make the median frame to do pixel replacement
            # not super efficient, but it works
            median_sf = mapwindow(median, sf.data, (median_size, median_size))

            # finally, assign the median values to the masked pixels
            sf.data[mask] .= median_sf[mask]

            reduced = AstroImage((sf .- matched_dark) ./ matched_flat, sf.header)

            reduced ./= reduced["COADDS"] # divide by coadds
            reduced .*= get_NIRC2_gain(reduced["DATE-OBS"]) # apply the gain
            
            reduced_filename = "reduced_$(lpad(reduced["FRAMENO"], 4, '0')).fits"
            reduced["RED-FN"] = reduced_filename

            reduced_frames[i] = reduced

        end

        @info "Reduced $(length(reduced_frames)) science frames"

        reduced_folder = joinpath(obslog["data_folder"], "reduced")
        make_and_clear(reduced_folder, "reduced_*.fits")

        reduced_filepaths = String[]
        for rf in reduced_frames
            reduced_filepath = joinpath(reduced_folder, rf["RED-FN"])
            push!(reduced_filepaths, rf["RED-FN"])
            @info "Saving reduced frame to $(reduced_filepath)"

            save(reduced_filepath, rf)
        end

        reduced_obslog = OrderedDict{String, Any}("data_folder" => obslog["data_folder"],
                                                  "date" => obslog["date"],
                                                  "reduced" => OrderedDict{String, Any}(
                                                  "reduced_sci" => reduced_filepaths))

        reduced_obslog_filepath = joinpath((obslog_folder, "$(obslog["date"])_reduced.toml"))

        @info "Writing to $(reduced_obslog_filepath)"
        toml_str = pretty_print_toml(reduced_obslog)
        open(reduced_obslog_filepath, "w") do io
            write(io, toml_str)
        end

        rejects_obslog_filepath = joinpath(obslog_folder, "$(obslog["date"])_rejects.toml")

        if isfile(rejects_obslog_filepath)
            @warn "Existing rejects file found! Skipping!"
        else
            rejects_obslog = OrderedDict{String, Any}("reduced" => OrderedDict{String, Any}("rejects" => String[]))
            @info "Writing to $(rejects_obslog_filepath)"
            toml_str = pretty_print_toml(rejects_obslog)
            open(rejects_obslog_filepath, "w") do io
                write(io, toml_str)
            end
        end

    end

end