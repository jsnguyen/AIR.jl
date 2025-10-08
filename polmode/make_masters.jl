using Glob
using Statistics
using AstroImages

using AIR

function make_darks(darks_frames; darks_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS"])

    if !isempty(darks_frames)
        master_darks, master_darks_masks = make_masters(darks_frames, darks_keylist, method=median)
    else
        master_darks = Dict{Any, AstroImage}()
        master_darks_masks = Dict{Any, BitMatrix}()
    end

    return collect(values(master_darks)), collect(values(master_darks_masks))
end

function make_master_darks_obslog(obslog; darks_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS"])
    darks_frames = load_frames(obslog, "darks")
    master_darks, master_darks_masks = make_darks(darks_frames; darks_keylist=darks_keylist)
    return master_darks, master_darks_masks
end

function make_flats(flat_frames, master_darks; flats_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")
    if !isempty(flat_frames)
        master_flats, master_flats_masks = make_masters(flat_frames, flats_keylist, method=median)

        for key in keys(master_flats)
            matched_dark = find_closest_dark(master_flats[key], master_darks)

            if matched_dark !== nothing
                @info "Subtracting dark from master flat"
                master_flats[key] .-= matched_dark
                master_flats[key]["FLATTYPE"] = flattype
            else
                @warn "No matching dark found for master flat"
                master_flats[key]["FLATTYPE"] = "$(flattype)+NODARK"
            end

            master_flats[key] ./= master_flats[key]["MAMEDIAN"]
        end

    else
        @warn "No regular flats found, skipping..."
        master_flats = Dict{Any, AstroImage}()
        master_flats_masks = Dict{Any, BitMatrix}()
    end
    return master_flats, master_flats_masks
end

function make_lamp_flats(lampon_frames, lampoff_frames, master_darks; flats_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"])
    if !isempty(lampon_frames)
        master_lampon, master_lampon_masks = make_masters(lampon_frames, flats_keylist, method=median)
    else
        @warn "No lampon frames found, skipping..."
        master_lampon = Dict{Any, AstroImage}()
        master_lampon_masks = Dict{Any, BitMatrix}()
    end

    if !isempty(lampoff_frames)
        master_lampoff, master_lampoff_masks = make_masters(lampoff_frames, flats_keylist, method=median)
    else
        @warn "No lampoff frames found, skipping..."
        master_lampoff = Dict{Any, AstroImage}()
        master_lampoff_masks = Dict{Any, BitMatrix}()
    end

    master_flat_lamp = Dict{Any, AstroImage}()
    master_flat_lamp_masks = Dict{Any, BitMatrix}()
    for key in keys(master_lampon)
        sub_frame = zeros(size(master_lampon[key]))

        master_flat_lamp[key] = master_lampon[key]
        
        # check if there is a corresponding lamp-off flat
        if haskey(master_lampoff, key)
            sub_frame = master_lampoff[key]
            master_flat_lamp[key]["FLATTYPE"] = "LAMP"
            master_flat_lamp_masks[key] = master_lampon_masks[key] .| master_lampoff_masks[key] # combine masks

        # if there is no corresponding lamp-off flat, use the dark frame
        else

            matched_dark = find_closest_dark(master_lampon[key], master_darks)
            if matched_dark === nothing
                @warn "No matching dark found for lamp-on flat"
                master_flat_lamp[key]["FLATTYPE"] = "LAMP+NODARK"
            else
                @warn "Using dark frame to subtract from lamp-on flat"
                sub_frame = matched_dark
                master_flat_lamp[key]["FLATTYPE"] = "LAMP+DARK"
            end
            master_flat_lamp_masks[key] = master_lampon_masks[key] # just use lampon mask if no lampoff frame is found

        end

        # subtract the lamp-off flat or dark frame from the lamp-on flat
        master_flat_lamp[key] .-= sub_frame

        # now get the median value
        # this should be bad pixel corrected already
        medval = median(master_flat_lamp[key])

        master_flat_lamp[key] ./= medval

    end

    return master_flat_lamp, master_flat_lamp_masks
end

function make_master_flats_obslog(obslog, master_darks, flats_keylist = ["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")

    #
    # regular flat frames
    #

    flat_frames = load_frames(obslog, "flats")
    master_flats, master_flats_masks = make_flats(flat_frames, master_darks; flats_keylist=flats_keylist)

    #
    # sky flat frames
    #

    sky_frames = load_frames(obslog, "flats_sky")
    master_flats_sky, master_flats_sky_masks = make_flats(sky_frames, master_darks; flats_keylist=flats_keylist, flattype="SKY")

    #
    # lamp flat frames
    #

    lampon_frames = load_frames(obslog, "flats_lampon")
    lampoff_frames = load_frames(obslog, "flats_lampoff")
    master_flat_lamp, master_flat_lamp_masks = make_lamp_flats(lampon_frames, lampoff_frames, master_darks; flats_keylist=flats_keylist)

    # ranking so that the best type of flat for a given key is always first on the list
    # when searching through flats, just pick the first one that matches
    master_flats = merge(master_flats, master_flats_sky, master_flat_lamp)
    master_flats = collect(values(master_flats)) # convert Dict to Vector
    ranking = ["SKY", "LAMP", "LAMP+DARK", "REGULAR", "SKY+NODARK", "LAMP+NODARK", "REGULAR+NODARK"]
    rank_dict = Dict(r => i for (i, r) in enumerate(ranking))
    sort!(master_flats, by = frame -> get(rank_dict, frame["FLATTYPE"], length(ranking)+1))

    master_flats_masks = merge(master_flats_masks, master_flats_sky_masks, master_flat_lamp_masks)
    master_flats_masks = collect(values(master_flats_masks)) # convert Dict to Vector

    return master_flats, master_flats_masks
end

function make_master_masks(master_darks_masks, master_flats_masks)
    # sort masks by size
    mask_stack = Dict{Any, Vector{BitMatrix}}()
    for mask in vcat(master_darks_masks, master_flats_masks)
        key = size(mask)

        if !haskey(mask_stack, key)
            mask_stack[key] = BitMatrix[]
        end

        push!(mask_stack[key], mask)

    end

    # combine masks by size
    master_masks = Dict{Any, Array{UInt8}}()
    for key in keys(mask_stack)
        mm = reduce(.|, mask_stack[key])
        mm = Array{UInt8}(mm) # convert BitMatrix to Array{UInt8} since FITS files don't support BitMatrix
        master_masks[key] = mm
    end

    return collect(values(master_masks))
end

@autolog begin

    obslog_folder = "polmode/obslogs"
    for obslog_filename in Glob.glob("*_obslog.toml", obslog_folder)

        @info "Loading obslog from" obslog_filename
        obslog = Obslog(obslog_filename)
        
        @info "Making darks..."
        master_darks, master_darks_masks = make_master_darks_obslog(obslog)
        @info "Writing master darks to" obslog.paths.darks_file
        save(obslog.paths.darks_file, master_darks...)

        @info "Making flats..."
        master_flats, master_flats_masks = make_master_flats_obslog(obslog, master_darks)
        @info "Writing master flats to" obslog.paths.flats_file
        save(obslog.paths.flats_file, master_flats...)

        @info "Making master masks..."
        master_masks = make_master_masks(master_darks_masks, master_flats_masks)
        @info "Writing masters masks to" obslog.paths.masks_file
        save(obslog.paths.masks_file, master_masks...)

    end

end
