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