using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LsqFit
using Plots
using ADI

using AIR

function imshow(img; kwargs...)
    xs, ys = axes(img)
    heatmap(xs, ys, transpose(img); aspect_ratio=1,
            xlim=extrema(xs), ylim=extrema(ys), kwargs...)
end;

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")

    cube = load(joinpath(sequences_folder, "aligned_frames.fits"), :)
    
    frames = cat([f.data for f in cube]...; dims=3)
    println("Loaded $(size(frames)[end]) frames from aligned_frames.fits")
    headers = [f.header for f in cube]

    alg = Classic() # 10 components
    angles = (x -> calculate_north_angle(x)[1]).(headers)
    reduced = alg(frames, angles)

    imshow(reduced)


end