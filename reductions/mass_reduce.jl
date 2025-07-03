using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics

using AIR

function load_master(master_filename)
    if isfile(master_filename)
        master = load(master_filename, :)
    else
        @warn "No master found..."
        master = AstroImage[]
    end
    return master
end

function load_masks(masks_filename)
    # get all the masks and load by size
    masks = Dict{Any, BitMatrix}()
    if isfile(masks_filename)
        ms = load(masks_filename, :)

        for m in ms 
            key = size(m)
            masks[key] = BitMatrix(m)
        end

    end
    return masks
end

autolog("$(@__FILE__).log") do

    obslog_folder = "reductions/obslogs"
    for obslog_filename in Glob.glob("*_obslog.toml", obslog_folder)
        @info "Loading obslog from" obslog_filename
        obslog = load_obslog(obslog_filename)
        
        darks_filename = joinpath(obslog["data_folder"], "darks.fits")
        flats_filename = joinpath(obslog["data_folder"], "flats.fits")
        masks_filename = joinpath(obslog["data_folder"], "master_mask.fits")

        master_darks = load_master(darks_filename)
        master_flats = load_master(flats_filename)
        masks = load_masks(masks_filename)

        sci_frames = load_frames(obslog, "sci")

        # add cosmic ray image cleaning here at some point?
        # first try didn't work properly...

        reduced_frames = AstroImage[]
        for sf in sci_frames

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
                mask = crop(NIRC2_bad_pixel_mask, size(sf))
            end

            if haskey(masks, size(sf))
                mask = mask .| masks[size(sf)] # also combine the extra mask if it exists
            end

            nan_mask = isnan.(sf.data) .| isinf.(sf.data) # create a mask for NaN and Inf values
            mask = mask .| nan_mask # combine the bad pixel mask with the NaN/Inf mask

            # make the median frame to do pixel replacement
            # not super efficient, but it works
            median_size = 5
            median_sf = mapwindow(median, sf.data, (median_size, median_size))

            # finally, assign the median values to the masked pixels
            sf.data[mask] .= median_sf[mask]

            reduced = AstroImage((sf .- matched_dark) ./ matched_flat, sf.header)

            reduced ./= reduced["COADDS"] # divide by coadds
            reduced .*= get_NIRC2_gain(reduced["DATE-OBS"]) # apply the gain

            push!(reduced_frames, reduced)

            reduced_filename = "reduced_$(lpad(reduced["FRAMENO"], 4, '0')).fits"
            reduced["RED-FN"] = reduced_filename

        end

        @info "Reduced $(length(reduced_frames)) science frames"

        reduced_folder = joinpath(obslog["data_folder"], "reduced")
        make_and_clear(reduced_folder, "reduced_*.fits")

        reduced_filepaths = String[]
        for rf in reduced_frames
            reduced_filepath = joinpath(reduced_folder, rf["RED-FN"])
            push!(reduced_filepaths, rf["RED-FN"])
            println("Saving reduced frame to $(reduced_filepath)")

            save(reduced_filepath, rf)
        end

        reduced_obslog = OrderedDict{String, Any}("data_folder" => obslog["data_folder"],
                                                  "subfolder" => "reduced",
                                                  "date" => obslog["date"],
                                                  "reduced" => reduced_filepaths)

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
            rejects_obslog = OrderedDict{String, Any}("rejects" => String[])
            @info "Writing to $(rejects_obslog_filepath)"
            toml_str = pretty_print_toml(rejects_obslog)
            open(rejects_obslog_filepath, "w") do io
                write(io, toml_str)
            end
        end

    end

end