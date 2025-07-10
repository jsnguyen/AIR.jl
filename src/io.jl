struct Obslog

    obslog::OrderedDict{String, Any}
    raw_folder::String
    reduced_folder::String
    plots_folder::String
    sequences_folder::String
    master_darks::Any
    master_flats::Any
    masks::Any
    sci::Vector{AstroImage}
    reduced_sci::Vector{AstroImage}
    sequences::Dict{String, Any}

    function Obslog(obslog_filename::String; rejects::Vector{String}=String[])
        obslog_dict = load_obslog(obslog_filename)
        
        darks_filename = joinpath(obslog_dict["data_folder"], "darks.fits")
        flats_filename = joinpath(obslog_dict["data_folder"], "flats.fits")
        masks_filename = joinpath(obslog_dict["data_folder"], "master_mask.fits")

        master_darks = load_master(darks_filename)
        master_flats = load_master(flats_filename)
        masks = load_masks(masks_filename)

        raw_folder = joinpath(obslog_dict["data_folder"], "raw")
        reduced_folder = joinpath(obslog_dict["data_folder"], "reduced")
        plots_folder = joinpath(obslog_dict["data_folder"], "plots")
        sequences_folder = joinpath(obslog_dict["data_folder"], "sequences")

        sci = load_frames(obslog_dict, "sci"; rejects=rejects)
        reduced_sci = load_frames(obslog_dict, "reduced_sci"; rejects=rejects)

        sequences = load_sequences(obslog_dict, rejects=rejects)

        new(obslog_dict, raw_folder, reduced_folder, plots_folder, sequences_folder, master_darks, master_flats, masks, sci, reduced_sci, sequences)
    end

end

Base.getindex(obs::Obslog, key) = obs.obslog[key]
Base.setindex!(obs::Obslog, value, key) = obs.obslog[key] = value
Base.keys(obs::Obslog) = keys(obs.obslog)
Base.values(obs::Obslog) = values(obs.obslog)
Base.haskey(obs::Obslog, key) = haskey(obs.obslog, key)
Base.length(obs::Obslog) = length(obs.obslog)
Base.iterate(obs::Obslog) = iterate(obs.obslog)
Base.iterate(obs::Obslog, state) = iterate(obs.obslog, state)
Base.get(obs::Obslog, key, default) = get(obs.obslog, key, default)

function load_obslog(obslog_filename::String)
    obslog = open(obslog_filename, "r") do file
        TOML.parse(file)
    end

    # unfold the nested dictionary to make things easier to access
    # the "raw" data is stored in a subfolder
    # the title of the .toml section refers to the subfolder where the files are located
    # note that if the folder name is the same as the key, this will break!

    for key in ["raw", "reduced"]
        if haskey(obslog, key)
            for k in keys(obslog[key])

                if k==key
                    error("Key '$(key)' cannot be the same as the folder name '$(k)'. Please rename the key or folder.")
                end

                obslog[k] = String[]
                for fn in obslog[key][k]
                    push!(obslog[k], joinpath(obslog["data_folder"], key, fn))
                end
            end
            delete!(obslog, key)
        end
    end

    return obslog
end

function load_frames(obslog, key; rejects::Vector{String}=String[])

    frames = AstroImage[]

    if !haskey(obslog, key)
        return frames
    end

    for fn in obslog[key]
        if basename(fn) in rejects
            continue
        end
        push!(frames, Float64.(load(fn)))
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

function load_sequences(sequence_obslog; rejects::Vector{String}=String[])
    sequences = Dict{String, Any}()
    for key in keys(sequence_obslog)
        if !(key in ["data_folder", "date"])
            sequences[key] = load_frames(sequence_obslog, key; rejects=rejects)
        end
    end
    return sequences
end

function load_rejects(rejects_path)
    rejects = open(rejects_path, "r") do file
        TOML.parse(file)["rejects"]
    end
    return rejects
end
