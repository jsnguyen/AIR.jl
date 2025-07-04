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
        sequences["$(name)_$(i)"] = seq[key]
    end
end

autolog("$(@__FILE__).log") do

    reduced_obslog_folder = "reductions/obslogs"
    reduced_obslog_path = joinpath(reduced_obslog_folder, "2002-06-16_reduced.toml")
    
    @info "Loading reduced_obslog from" reduced_obslog_path
    reduced_obslog = load_obslog(reduced_obslog_path)
    rejects_obslog_path = joinpath(reduced_obslog["data_folder"], reduced_obslog["subfolder"], "rejects.toml")
    rejects_obslog = load_obslog(rejects_obslog_path)
    reduced = load_frames(reduced_obslog, "reduced", rejects=rejects_obslog["rejects"])

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

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog["data_folder"],
                                       "subfolder" => "reduced",
                                       "date" => reduced_obslog["date"])

    add_to_sequence!(sequences, hbc630_seq, "hbc630")
    add_to_sequence!(sequences, hbc650_seq, "hbc650")
    add_to_sequence!(sequences, as209_seq, "as209")

    sequence_filepath = joinpath(reduced_obslog_folder, "$(reduced_obslog["date"])_sequences.toml")
    toml_str = pretty_print_toml(sequences)
    open(sequence_filepath, "w") do io
        write(io, toml_str)
    end

end