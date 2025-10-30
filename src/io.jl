struct ObslogPaths
    date::String

    data_folder::String

    raw_folder::String
    reduced_folder::String
    plots_folder::String
    sequences_folder::String

    reduced_file::String
    rejects_file::String
    sequences_file::String
    table_file::String

    darks_file::String
    flats_file::String
    skies_file::String
    masks_file::String

    function ObslogPaths(date::String, data_folder::String)

        raw_folder = joinpath(data_folder, "raw")
        reduced_folder = joinpath(data_folder, "reduced")
        plots_folder = joinpath(data_folder, "plots")
        sequences_folder = joinpath(data_folder, "sequences")

        reduced_file = joinpath(data_folder, "$(date)_reduced.toml")
        rejects_file = joinpath(data_folder, "$(date)_rejects.toml")
        sequences_file = joinpath(data_folder, "$(date)_sequences.toml")
        table_file = joinpath(data_folder, "$(date)_reduced_frames_table.txt")

        darks_file = joinpath(data_folder, "darks.fits")
        flats_file = joinpath(data_folder, "flats.fits")
        skies_file = joinpath(data_folder, "skies.fits")
        masks_file = joinpath(data_folder, "master_mask.fits")

        new(date, data_folder,
            raw_folder, reduced_folder, plots_folder, sequences_folder,
            reduced_file, rejects_file, sequences_file, table_file,
            darks_file, flats_file, skies_file, masks_file)
    end

    function ObslogPaths(obslog_dict::OrderedDict{String, Any}, date)
        data_folder = obslog_dict["data_folder"]
        delete!(obslog_dict, "data_folder")

        raw_folder = joinpath(data_folder, "raw")
        reduced_folder = joinpath(data_folder, "reduced")
        plots_folder = joinpath(data_folder, "plots")
        sequences_folder = joinpath(data_folder, "sequences")

        reduced_file = joinpath(data_folder, "$(date)_reduced.toml")
        rejects_file = joinpath(data_folder, "$(date)_rejects.toml")
        sequences_file = joinpath(data_folder, "$(date)_sequences.toml")
        table_file = joinpath(data_folder, "$(date)_reduced_frames_table.txt")

        darks_file = joinpath(data_folder, "darks.fits")
        flats_file = joinpath(data_folder, "flats.fits")
        skies_file = joinpath(data_folder, "skies.fits")
        masks_file = joinpath(data_folder, "master_mask.fits")

        new(date, data_folder,
            raw_folder, reduced_folder, plots_folder, sequences_folder,
            reduced_file, rejects_file, sequences_file, table_file,
            darks_file, flats_file, skies_file, masks_file)
    end

end

struct Obslog

    obslog::OrderedDict{String, Any}
    paths::ObslogPaths

    date::String

    master_darks::Any
    master_flats::Any
    master_skies::Any
    masks::Any
    sci::Vector{AstroImage}
    reduced_sci::Vector{AstroImage}
    sequences::Dict{String, Any}

    rejects::Vector{String}

    is_loaded::Bool

    function Obslog(obslog_path::String)
        obslog_dict = load_obslog_dict_from_file(obslog_path)
        obslog = _make_obslog(obslog_dict)
        return obslog
    end

    function Obslog(obslog_dict::OrderedDict{String, Any})
        obslog = _make_obslog(obslog_dict)
        return obslog
    end

    function _make_obslog(obslog_dict::OrderedDict{String, Any})

        unfold_obslog_dict(obslog_dict)

        # shortcuts
        date = obslog_dict["date"]
        delete!(obslog_dict, "date")

        paths = ObslogPaths(obslog_dict, date)

        loaded_files = _load_files(obslog_dict, paths)

        return new(obslog_dict, paths, date, loaded_files...)

    end

    function _load_files(obslog_dict::OrderedDict{String, Any}, paths::ObslogPaths)
        rejects = String[]
        if isfile(paths.rejects_file)
            @info "Rejects file found: $(paths.rejects_file)"
            rejects = load_rejects(paths.rejects_file)
        end

        master_darks = load_master(paths.darks_file)
        master_flats = load_master(paths.flats_file)
        master_skies = load_master(paths.skies_file)
        masks = load_masks(paths.masks_file)

        sci = load_frames(obslog_dict, "sci"; rejects=rejects)
        reduced_sci = load_frames(obslog_dict, "reduced_sci"; rejects=rejects)

        protected_keys = ["obslog_path"]
        sequences = Dict{String, Any}()
        if occursin("sequences", obslog_dict["obslog_path"])
            for key in keys(obslog_dict)
                if !(key in protected_keys)
                    sequences[key] = load_frames(obslog_dict, key; rejects=rejects)
                end
            end
        end

        is_loaded = true

        return master_darks, master_flats, master_skies, masks, sci, reduced_sci, sequences, rejects, is_loaded

    end

    function _unload_files!(obslog::Obslog)
        obslog.master_darks = Dict{Any, AstroImage}()
        obslog.master_flats = Dict{Any, AstroImage}()
        obslog.master_skies = Dict{Any, AstroImage}()

        obslog.masks = Dict{Any, BitMatrix}()
        obslog.sci = Dict{Any, AstroImage}()
        obslog.reduced_sci = Dict{Any, AstroImage}()

        obslog.sequences = Dict{String, Any}()
        obslog.rejects = String[]

        obslog.is_loaded = false
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

function unfold_obslog_dict(obslog_dict::OrderedDict{String, Any})
    # unfold the nested dictionary to make things easier to access
    # the "raw" data is stored in a subfolder
    # the title of the .toml section refers to the subfolder where the files are located
    # note that if the folder name is the same as the key, this will break!

    for key in ["raw", "reduced"]
        if haskey(obslog_dict, key)
            for k in keys(obslog_dict[key])

                if k==key
                    error("Key '$(key)' cannot be the same as the folder name '$(k)'. Please rename the key or folder.")
                end

                obslog_dict[k] = String[]
                for fn in obslog_dict[key][k]
                    push!(obslog_dict[k], joinpath(obslog_dict["data_folder"], key, fn))
                end
            end
            delete!(obslog_dict, key)
        end
    end

    return obslog_dict
end

function load_obslog_dict_from_file(obslog_path::String)
    obslog_dict = open(obslog_path, "r") do file
        TOML.parse(file)
    end
    obslog_dict["obslog_path"] = obslog_path
    obslog_dict = OrderedDict(pairs(obslog_dict))

    return obslog_dict
end

function load_frames(obslog::Union{Obslog, OrderedDict{String, Any}}, key; rejects::Vector{String}=String[])

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

function load_frames(frame_paths::Vector{String}; rejects::Vector{String}=String[])

    frames = AstroImage[]

    for fn in frame_paths
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