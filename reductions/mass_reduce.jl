using Glob
using AstroImages
using ImageFiltering
using Statistics

using AIR


autolog("$(@__FILE__).log") do

    for obslog_filename in Glob.glob("reductions/obslogs/*.toml")
        @info "Loading obslog from" obslog_filename
        obslog = load_obslog(obslog_filename)
        
        darks_filename = joinpath(obslog["data_folder"], "darks.fits")
        flats_filename = joinpath(obslog["data_folder"], "flats.fits")
        masks_filename = joinpath(obslog["data_folder"], "master_mask.fits")

        if isfile(darks_filename)
            master_darks = load(darks_filename, :)
        else
            @warn "No master darks found..."
            master_darks = AstroImage[]
        end

        if isfile(flats_filename)
            master_flats = load(flats_filename, :)
        else
            @warn "No master flats found..."
            master_flats = AstroImage[]
        end

        masks = Dict{Any, BitMatrix}()
        if isfile(masks_filename)
            raw_masks = load(masks_filename, :)

            for rm in raw_masks
                key = size(rm)
                masks[key] = BitMatrix(rm)
            end

        end

        sci_frames = load_frames(obslog, "sci")

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

            bad_pixel_mask = copy(NIRC2_bad_pixel_mask)
            if size(bad_pixel_mask) != size(sf)
                bad_pixel_mask = crop(NIRC2_bad_pixel_mask, size(sf))
            end

            mask = bad_pixel_mask .| masks[size(sf)] # also combine the extra mask if it exists

            # make the median frame to do pixel replacement
            # not super efficient, but it works
            median_size = 5
            median_sf = mapwindow(median, sf.data, (median_size, median_size))

            # finally, assign the median values to the masked pixels
            sf.data[mask] .= median_sf[mask]

            reduced = AstroImage((sf .- matched_dark) ./ matched_flat, sf.header)
            push!(reduced_frames, reduced)

        end

        @info "Reduced $(length(reduced_frames)) science frames"

        reduced_folder = joinpath(obslog["data_folder"], "reduced")
        if !isdir(reduced_folder)
            mkpath(reduced_folder)
        end

        for (i,rf) in enumerate(reduced_frames) 
            reduced_filename = joinpath(reduced_folder, "reduced_$(lpad(string(i), 4, '0')).fits")
            println("Saving reduced frame to $(reduced_filename)")

            save(reduced_filename, rf)
        end

    end

end
