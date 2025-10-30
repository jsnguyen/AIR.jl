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
        @info "adding sequence" name i key length(seq[key])
        sequences["reduced"]["$(name)_$(i)"] = seq[key]
    end
end

@stage function generate_sequence_epoch(reduced_frames, indices_between)

    # got these from looking at the framelist and checking the frames
    obj_frames = get_between(reduced_frames, indices_between)
    @info "Found $(length(obj_frames)) object frames between $(indices_between)..."
    keylist = ["ITIME", "FILTER"]
    obj_seq = match_keys(obj_frames, keylist)
    sequences = collect(values(obj_seq))

    return sequences

end