function load_obslog(obslog_filename::String)
    obslog = open(obslog_filename, "r") do file
        TOML.parse(file)
    end
    return obslog
end

function load_frames(obslog, key; rejects::Vector{String}=String[])

    frames = AstroImage[]

    if !haskey(obslog, key)
        return frames
    end

    for fn in obslog[key]
        if fn in rejects
            continue
        end
            push!(frames, Float64.(load(joinpath(obslog["data_folder"], obslog["subfolder"], fn))))
    end

    return frames
end

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