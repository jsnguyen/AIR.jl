using OrderedCollections
using AstroImages
using Statistics

using AIR

get_between(frames, between) = filter(frame -> between[1] <= frame["FRAMENO"] <= between[2], frames)

function sequence_dict(frames, keylist)
    matched_dict = match_keys(frames, keylist)

    str_dict = Dict{Any,Vector{String}}()
    for key in keys(matched_dict)
        str_dict[key] = []
        for f in matched_dict[key]
            if !isnothing(f)
                push!(str_dict[key], f["RED-FN"])
            end
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

function empty_sequence(data_folder, date, seq_type)
    OrderedDict{String,Any}("data_folder" => data_folder,
                             "date" => date,
                             seq_type => OrderedDict{String,Any}())
end

# 2002-06-16, Epoch 1
function generate_sequence_epoch_1()

    obslog_folder = "pipeline/obslogs"
    reduced_obslog_path = joinpath(obslog_folder, "2002-06-16_reduced.toml")
    rejects_obslog_path = joinpath(obslog_folder, "2002-06-16_rejects.toml")
    rejects = load_rejects(rejects_obslog_path)

    @info "Loading reduced_obslog from" reduced_obslog_path
    reduced_obslog = Obslog(reduced_obslog_path, rejects=rejects)
    reduced = reduced_obslog.reduced_sci

    # got these from looking at the framelist and checking the frames
    hbc630_between = (412, 440)
    hbc650_between = (441, 465)
    as209_between = (466, 490)

    hbc630_frames = get_between(reduced, hbc630_between)
    hbc650_frames = get_between(reduced, hbc650_between)
    as209_frames = get_between(reduced, as209_between)
    
    keylist = ["ITIME", "FILTER"]
    hbc630_seq = sequence_dict(hbc630_frames, keylist)
    hbc650_seq = sequence_dict(hbc650_frames, keylist)
    as209_seq = sequence_dict(as209_frames, keylist)

    sequences = empty_sequence(reduced_obslog["data_folder"], reduced_obslog["date"], "reduced")

    add_to_sequence!(sequences, hbc630_seq, "hbc630")
    add_to_sequence!(sequences, hbc650_seq, "hbc650")
    add_to_sequence!(sequences, as209_seq, "as209")

    sequence_filepath = joinpath(reduced_obslog_folder, "$(reduced_obslog["date"])_sequences.toml")
    toml_str = pretty_print_toml(sequences)
    open(sequence_filepath, "w") do io
        write(io, toml_str)
    end

    return sequences

end

@autolog begin

    generate_sequence_epoch_1()

end