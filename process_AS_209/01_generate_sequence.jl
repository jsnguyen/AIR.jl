using OrderedCollections
using AstroImages
using Statistics

using AIR

get_between(frames, between) = filter(frame -> between[1] <= frame["FRAMENO"] <= between[2], frames)

function sequence_dict(frames, keylist; ignore=Dict{String,Any}())
    matched_dict = match_keys(frames, keylist)

    str_dict = Dict{Any,Vector{AstroImage}}()
    for key in keys(matched_dict)
        str_dict[key] = []
        for f in matched_dict[key]

            ignore_flag = false
            for k in keys(ignore)
                if f[k] == ignore[k]
                    @info "Ignoring frame $(f["RED-FN"]) due to ignore condition: $k = $(ignore[k])"
                    ignore_flag = true
                    break
                end
            end
            
            if !isnothing(f) && !ignore_flag
                push!(str_dict[key], f)
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
        sequences["$(name)_$(i)"] = seq[key]
    end
end

@stage function generate_sequence_epoch(targets, reduced_frames)

    sequences = Dict{String, Any}()
    for (targ, inds) in targets

        obj_frames = get_between(reduced_frames, inds)
        @info "Found $(length(obj_frames)) object frames between $(inds)..."
        keylist = ["ITIME", "FILTER"]
        obj_seq = sequence_dict(obj_frames, keylist; ignore=Dict("CAMNAME" => "wide"))

        @info "Generated sequences for target: $targ"
        for (i,key) in enumerate(keys(obj_seq))
            @info "Sequence key $(i): $(key) has $(length(obj_seq[key])) frames"
        end

        for (i,key) in enumerate(keys(obj_seq))
            @info targ=targ inds=inds i=i key=key length=length(obj_seq[key])
        end

        seq = collect(values(obj_seq))
        add_to_sequence!(sequences, seq, targ)
    end

    @info "Generated sequences: keys=$(sort(collect(keys(sequences))))"

    return sequences

end