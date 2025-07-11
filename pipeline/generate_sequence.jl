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

# 2002-06-16, Epoch 1
function generate_sequence_epoch_1()

    reduced_obslog = Obslog("pipeline/obslogs/2002-06-16_reduced.toml")
    @info "Loading reduced_obslog from" reduced_obslog.paths.obslog_file

    # got these from looking at the framelist and checking the frames
    hbc630_between = (412, 440)
    hbc650_between = (441, 465)
    as209_between = (466, 490)

    hbc630_frames = get_between(reduced_obslog.reduced_sci, hbc630_between)
    hbc650_frames = get_between(reduced_obslog.reduced_sci, hbc650_between)
    as209_frames = get_between(reduced_obslog.reduced_sci, as209_between)
    
    keylist = ["ITIME", "FILTER"]
    hbc630_seq = sequence_dict(hbc630_frames, keylist)
    hbc650_seq = sequence_dict(hbc650_frames, keylist)
    as209_seq = sequence_dict(as209_frames, keylist)

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog.paths.data_folder,
                                        "date" => reduced_obslog.date,
                                        "reduced" => OrderedDict{String,Any}())

    add_to_sequence!(sequences, hbc630_seq, "hbc630")
    add_to_sequence!(sequences, hbc650_seq, "hbc650")
    add_to_sequence!(sequences, as209_seq, "as209")

    write_toml(reduced_obslog.paths.sequences_file, sequences)

    return sequences

end

# 2002-08-02, Epoch 2
function generate_sequence_epoch_2()

    reduced_obslog = Obslog("pipeline/obslogs/2002-08-02_reduced.toml")
    @info "Loading reduced_obslog from" reduced_obslog.paths.obslog_file

    # got these from looking at the framelist and checking the frames
    # skip first 2 for being in wrong filter
    as209_between = (78, 109)

    as209_frames = get_between(reduced_obslog.reduced_sci, as209_between)
    
    keylist = ["ITIME", "FILTER"]
    as209_seq = sequence_dict(as209_frames, keylist; ignore=Dict("CAMNAME" => "wide"))

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog.paths.data_folder,
                                        "date" => reduced_obslog.date,
                                        "reduced" => OrderedDict{String,Any}())

    add_to_sequence!(sequences, as209_seq, "as209")

    write_toml(reduced_obslog.paths.sequences_file, sequences)

    return sequences

end

# 2002-08-02, Epoch 3
function generate_sequence_epoch_3()

    reduced_obslog = Obslog("pipeline/obslogs/2002-08-21_reduced.toml")
    @info "Loading reduced_obslog from" reduced_obslog.paths.obslog_file

    # got these from looking at the framelist and checking the frames
    # skip first 2 for being in wrong filter
    as209_between = (84, 117)
    tyc2307_between = (192, 219)

    as209_frames = get_between(reduced_obslog.reduced_sci, as209_between)
    tyc2307_between_frames = get_between(reduced_obslog.reduced_sci, tyc2307_between)
    
    keylist = ["ITIME", "FILTER"]
    as209_seq = sequence_dict(as209_frames, keylist; ignore=Dict("CAMNAME" => "wide"))
    tyc2307_seq = sequence_dict(tyc2307_between_frames, keylist; ignore=Dict("CAMNAME" => "wide"))

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog.paths.data_folder,
                                        "date" => reduced_obslog.date,
                                        "reduced" => OrderedDict{String,Any}())

    add_to_sequence!(sequences, as209_seq, "as209")
    add_to_sequence!(sequences, tyc2307_seq, "tyc2307")

    write_toml(reduced_obslog.paths.sequences_file, sequences)

    return sequences

end

# 2005-07-27, Epoch 4
function generate_sequence_epoch_4()

    reduced_obslog = Obslog("pipeline/obslogs/2005-07-27_reduced.toml")
    @info "Loading reduced_obslog from" reduced_obslog.paths.obslog_file

    # got these from looking at the framelist and checking the frames
    # skip first 2 for being in wrong filter
    as209_between = (181, 205)
    t222007 = (373, 398)

    as209_frames = get_between(reduced_obslog.reduced_sci, as209_between)
    t222007_frames = get_between(reduced_obslog.reduced_sci, t222007)

    keylist = ["ITIME", "FILTER"]
    as209_seq = sequence_dict(as209_frames, keylist; ignore=Dict("CAMNAME" => "wide"))
    t222007_seq = sequence_dict(t222007_frames, keylist; ignore=Dict("CAMNAME" => "wide"))

    sequences = OrderedDict{String,Any}("data_folder" => reduced_obslog.paths.data_folder,
                                        "date" => reduced_obslog.date,
                                        "reduced" => OrderedDict{String,Any}())

    add_to_sequence!(sequences, as209_seq, "as209")
    add_to_sequence!(sequences, t222007_seq, "t222007")

    write_toml(reduced_obslog.paths.sequences_file, sequences)

    return sequences

end


@autolog begin

    generate_sequence_epoch_1()
    generate_sequence_epoch_2()
    generate_sequence_epoch_3()
    generate_sequence_epoch_4()

end