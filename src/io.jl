struct ObslogPaths
    obslog_folder::String
    data_folder::String

    raw_folder::String
    reduced_folder::String
    plots_folder::String
    sequences_folder::String

    obslog_file::String
    reduced_file::String
    rejects_file::String
    sequences_file::String
    table_file::String

    darks_file::String
    flats_file::String
    masks_file::String

    function ObslogPaths(obslog_dict::Dict{String, Any}, date)
        obslog_folder = dirname(obslog_dict["obslog_file"])
        data_folder = obslog_dict["data_folder"]
        delete!(obslog_dict, "data_folder")

        raw_folder = joinpath(data_folder, "raw")
        reduced_folder = joinpath(data_folder, "reduced")
        plots_folder = joinpath(data_folder, "plots")
        sequences_folder = joinpath(data_folder, "sequences")

        obslog_file = obslog_dict["obslog_file"]
        reduced_file = joinpath(obslog_folder, "$(date)_reduced.toml")
        rejects_file = joinpath(obslog_folder, "$(date)_rejects.toml")
        sequences_file = joinpath(obslog_folder, "$(date)_sequences.toml")
        table_file = joinpath(obslog_folder, "$(date)_reduced_frames_table.txt")

        darks_file = joinpath(data_folder, "darks.fits")
        flats_file = joinpath(data_folder, "flats.fits")
        masks_file = joinpath(data_folder, "master_mask.fits")

        new(obslog_folder, data_folder,
            raw_folder, reduced_folder, plots_folder, sequences_folder,
            obslog_file, reduced_file, rejects_file, sequences_file, table_file,
            darks_file, flats_file, masks_file)
    end

end

struct Obslog

    obslog::OrderedDict{String, Any}
    paths::ObslogPaths

    date::String

    master_darks::Any
    master_flats::Any
    masks::Any
    sci::Vector{AstroImage}
    reduced_sci::Vector{AstroImage}
    sequences::Dict{String, Any}

    function Obslog(obslog_file::String)
        obslog_dict = load_obslog(obslog_file)

        # shortcuts
        date = obslog_dict["date"]
        delete!(obslog_dict, "date")

        paths = ObslogPaths(obslog_dict, date)

        rejects = String[]
        if isfile(paths.rejects_file)
            @info "Rejects file found: $(paths.rejects_file)"
            rejects = load_rejects(paths.rejects_file)
        end
        
        master_darks = load_master(paths.darks_file)
        master_flats = load_master(paths.flats_file)
        masks = load_masks(paths.masks_file)

        sci = load_frames(obslog_dict, "sci"; rejects=rejects)
        reduced_sci = load_frames(obslog_dict, "reduced_sci"; rejects=rejects)

        protected_keys = ["obslog_file"]
        sequences = Dict{String, Any}()
        if occursin("sequences", obslog_file)
            for key in keys(obslog_dict)
                if !(key in protected_keys)
                    sequences[key] = load_frames(obslog_dict, key; rejects=rejects)
                end
            end
        end

        new(obslog_dict, paths, date, master_darks, master_flats, masks, sci, reduced_sci, sequences)
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

function load_obslog(obslog_file::String)
    obslog = open(obslog_file, "r") do file
        TOML.parse(file)
    end
    obslog["obslog_file"] = obslog_file

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

function load_rejects(rejects_file)
    rejects = open(rejects_file, "r") do file
        TOML.parse(file)["rejects"]
    end
    return rejects
end

function write_toml(filename, data_dict)
    @info "Writing sequences to" filename
    toml_str = pretty_print_toml(data_dict)
    open(filename, "w") do io
        write(io, toml_str)
    end
end
