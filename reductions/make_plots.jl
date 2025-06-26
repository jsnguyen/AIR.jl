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

        @info "Loading obslog from" obslog_filename
        obslog = load_obslog(obslog_filename)
        reduced = load_frames(obslog, "reduced")

        #
        # Plot
        #

        plotting_folder = joinpath(obslog["data_folder"], "plots")
        make_and_clear(plotting_folder, "frames_*.png")

        @showprogress for frame in reduced
            save(joinpath(plotting_folder, "frames_$(lpad(string(frame["FRAMENO"]), 4, '0')).png"), imview(frame, cmap=:matter))
        end

    end



end