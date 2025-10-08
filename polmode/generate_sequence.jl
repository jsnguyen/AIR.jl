using OrderedCollections
using AstroImages
using Statistics

using AIR

get_between(frames, between) = filter(frame -> between[1] <= frame["FRAMENO"] <= between[2], frames)

function sequence_dict(frames, keylist; ignore=Dict{String,Any}())
    matched_dict = match_keys(frames, keylist)

    str_dict = Dict{Any,Vector{String}}()
    for key in keys(matched_dict)
        str_dict[key] = []
        for f in matched_dict[key]

            ignore_flag = false
            for k in keys(ignore)
                @info "Checking ignore condition for key: $k with value $(ignore[k]) $(f[k])"
                if f[k] == ignore[k]
                    @info "Ignoring frame $(f["RED-FN"]) due to ignore condition: $k = $(ignore[k])"
                    ignore_flag = true
                    break
                end
            end
            
            if !isnothing(f) && !ignore_flag
                push!(str_dict[key], f["RED-FN"])
            end
        end

        # might be empty if no frames matched or frames got ignored
        if str_dict[key] == []
            delete!(str_dict, key)
        end

    end
    return str_dict
end


function add_to_sequence!(sequences, seq, name)
    for (i,key) in enumerate(keys(seq))
        @info "adding sequence" name i key length(seq[key])
        sequences["reduced"]["$(name)_$(i)"] = seq[key]
    end
end

function generate_sequence_epoch()

    reduced_obslog = Obslog("live_ingestion/obslogs/2025-10-07_reduced.toml")
    @info "Loading reduced_obslog from" reduced_obslog.paths.obslog_file

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog.paths.data_folder,
                                        "date" => reduced_obslog.date,
                                        "reduced" => OrderedDict{String,Any}())

    # got these from looking at the framelist and checking the frames
    obj_between = (247, 259)
    obj_frames = get_between(reduced_obslog.reduced_sci, obj_between)
    keylist = ["ITIME", "FILTER"]
    obj_seq = sequence_dict(obj_frames, keylist)
    add_to_sequence!(sequences, obj_seq, "HD215806")

    # got these from looking at the framelist and checking the frames
    obj_between = (276, 336)
    obj_frames = get_between(reduced_obslog.reduced_sci, obj_between)
    keylist = ["ITIME", "FILTER"]
    obj_seq = sequence_dict(obj_frames, keylist)
    add_to_sequence!(sequences, obj_seq, "BD45598")

    write_toml(reduced_obslog.paths.sequences_file, sequences)

    return sequences

end

@autolog begin

    generate_sequence_epoch()

end