using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics
using ProgressMeter

using AIR

autolog("$(@__FILE__).log") do

    reduced_obslog_folder = "reductions/obslogs"
    reduced_obslog_path = joinpath(reduced_obslog_folder, "2002-06-16_reduced.toml")

    @info "Loading reduced_obslog from" reduced_obslog_path
    reduced_obslog = load_obslog(reduced_obslog_path)
    reduced = load_frames(reduced_obslog, "reduced")

    unsat = String[]
    kp_seq = String[]
    j_seq = String[]
    for frame in reduced
        if occursin("PK50", frame["FILTER"])
            push!(unsat, frame["RED-FN"])
        elseif occursin("J", frame["FILTER"])
            push!(j_seq, frame["RED-FN"])
        elseif occursin("Kp", frame["FILTER"])
            push!(kp_seq, frame["RED-FN"])
        end
    end

    sequence = OrderedDict{String,Any}("data_folder" => reduced_obslog["data_folder"],
                                       "subfolder" => "reduced",
                                       "date" => reduced_obslog["date"],
                                       "unsat" => unsat,
                                       "j_seq" => j_seq,
                                       "kp_seq" => kp_seq)

    sequence_filepath = joinpath(reduced_obslog_folder, "$(reduced_obslog["date"])_sequences.toml")
    toml_str = pretty_print_toml(sequence)
    open(sequence_filepath, "w") do io
        write(io, toml_str)
    end

end