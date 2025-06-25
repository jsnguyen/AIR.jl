using OrderedCollections
using Glob
using AstroImages
using ImageFiltering
using Statistics
using ProgressMeter

using AIR


autolog("$(@__FILE__).log") do

    obslog_folder = "reductions/obslogs"
    for obslog_filename in Glob.glob(joinpath(obslog_folder,"*_reduced.toml"))

        #obslog_path = joinpath(obslog_folder, "2002-08-21_AS_209.toml")

        @info "Loading obslog from" obslog_filename
        obslog = load_obslog(obslog_filename)
        reduced_folder = joinpath(obslog["data_folder"], "reduced")

        frames = [load(fn) for fn in Glob.glob("*.fits", reduced_folder)]

        #
        # Plot
        #

        plotting_folder = joinpath(obslog["data_folder"], "plots")
        make_and_clear(plotting_folder, "frames_*.png")

        @showprogress for f in frames
            save(joinpath(plotting_folder, "frames_$(lpad(string(f["FRAMENO"]), 4, '0')).png"), imview(f, cmap=:matter))
        end

    end



end