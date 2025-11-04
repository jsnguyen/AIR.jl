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

function make_master_darks(darks_frames; darks_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS"])
    master_darks, master_darks_masks = make_darks(darks_frames; darks_keylist=darks_keylist)
    return master_darks, master_darks_masks
end

function make_flats(flat_frames, master_darks; flats_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")
    if !isempty(flat_frames)
        master_flats, master_flats_masks = make_masters(flat_frames, flats_keylist, method=median)

        for key in keys(master_flats)
            ind, matched_dark = find_closest_dark(master_flats[key], master_darks)

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

function make_skies(sky_frames, master_darks; flats_keylist=["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")
    if !isempty(sky_frames)
        master_skies, master_skies_masks = make_masters(sky_frames, flats_keylist, method=median)

        for key in keys(master_skies)
            ind, matched_dark = find_closest_dark(master_skies[key], master_darks)

            if matched_dark !== nothing
                @info "Subtracting dark from master flat"
                master_skies[key] .-= matched_dark
                master_skies[key]["FLATTYPE"] = flattype
            else
                @warn "No matching dark found for master flat"
                master_skies[key]["FLATTYPE"] = "$(flattype)+NODARK"
            end

        end

    else
        @warn "No regular flats found, skipping..."
        master_skies = Dict{Any, AstroImage}()
        master_skies_masks = Dict{Any, BitMatrix}()
    end

    master_skies = collect(values(master_skies))
    master_skies_masks = collect(values(master_skies_masks))
    return master_skies, master_skies_masks
end

function make_master_skies(sky_frames, master_darks, flats_keylist = ["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")
    master_skies, master_skies_masks = make_skies(sky_frames, master_darks; flats_keylist=flats_keylist, flattype=flattype)
    return master_skies, master_skies_masks
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

            ind, matched_dark = find_closest_dark(master_lampon[key], master_darks)
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

function make_master_flats(flat_frames, sky_frames, lampon_frames, lampoff_frames, master_darks, flats_keylist = ["NAXIS1", "NAXIS2", "ITIME", "COADDS", "FILTER"], flattype="REGULAR")

    #
    # regular flat frames
    #

    master_flats, master_flats_masks = make_flats(flat_frames, master_darks; flats_keylist=flats_keylist)

    #
    # sky flat frames
    #

    master_flats_sky, master_flats_sky_masks = make_flats(sky_frames, master_darks; flats_keylist=flats_keylist, flattype="SKY")

    #
    # lamp flat frames
    #

    master_flat_lamp, master_flat_lamp_masks = make_lamp_flats(lampon_frames, lampoff_frames, master_darks; flats_keylist=flats_keylist)

    # ranking so that the best type of flat for a given key is always first on the list
    # when searching through flats, just pick the first one that matches
    master_flats = merge(master_flats, master_flats_sky, master_flat_lamp)
    master_flats = collect(values(master_flats)) # convert Dict to Vector
    ranking = ["SKY", "LAMP", "LAMP+DARK", "REGULAR", "SKY+NODARK", "LAMP+NODARK", "REGULAR+NODARK"]
    rank_dict = Dict(r => i for (i, r) in enumerate(ranking))

    # more frames = better, so sort by the number of frames too, negative for descending order
    sort!(master_flats, by = frame -> (get(rank_dict, frame["FLATTYPE"], length(ranking)+1), -frame["NFRAMES"]))

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
    master_masks = Dict{Any, BitMatrix}()
    for key in keys(mask_stack)
        mm = reduce(.|, mask_stack[key])
        master_masks[key] = mm
    end

    return master_masks
end

@stage function make_all_masters(sci_filepaths, flats_filepaths, flats_sky_filepaths, flats_lampon_filepaths, flats_lampoff_filepaths, darks_filepaths)

    if isfile(context["paths"].rejects_file)
        rejects = load_rejects(context["paths"].rejects_file)
    else
        rejects = String[]
    end
    #sci = load.(sci_filepaths)
    flats = load_frames(flats_filepaths; rejects=rejects)
    flats_sky = load_frames(flats_sky_filepaths; rejects=rejects)
    flats_lampon = load_frames(flats_lampon_filepaths; rejects=rejects)
    flats_lampoff = load_frames(flats_lampoff_filepaths; rejects=rejects)
    darks = load_frames(darks_filepaths; rejects=rejects)

    paths = context["paths"]

    @info "Making darks..."
    master_darks, master_darks_masks = make_master_darks(darks)
    @info "Writing master darks to" paths.darks_file
    save(paths.darks_file, master_darks...)
    @context_store context paths.darks_file

    @info "Making flats..."
    master_flats, master_flats_masks = make_master_flats(flats, flats_sky, flats_lampon, flats_lampoff, master_darks)
    @info "Writing master flats to" paths.flats_file
    save(paths.flats_file, master_flats...)
    @context_store context paths.flats_file

    @info "Making skies..."
    master_skies, master_skies_masks = make_master_skies(flats_sky, master_darks)
    @info typeof(master_skies)
    @info "Writing master skies to" paths.skies_file
    save(paths.skies_file, master_skies...)
    @context_store context paths.skies_file

    @info "Making master masks..."
    master_masks = make_master_masks(master_darks_masks, master_flats_masks)
    @info "Writing masters masks to" paths.masks_file
    save(paths.masks_file, [Array{UInt8}(mm) for mm in collect(values(master_masks))]...)
    @context_store context paths.masks_file

    return master_darks, master_flats, master_skies, master_masks


end