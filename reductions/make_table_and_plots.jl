using Glob
using Printf
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

        fields = ["RED-FN", "OBJECT", "TARGNAME", "RA", "DEC", "CAMNAME", "DATE-OBS", "UTC", "ITIME", "COADDS", "FILTER", "NAXIS1", "NAXIS2"]
        field_lengths = [18, 16, 16, 12, 12, 8, 10, 11, 8, 8, 20, 6, 6]

        table_filename = joinpath(obslog_folder, "$(obslog["date"])_frames_table.txt")
        @info "Writing frames table to" table_filename

        open(table_filename, "w") do io
            # print header
            for (i, (l,f)) in enumerate(zip(field_lengths,fields))
                @printf(io, "%-*s ", l, f)
                if i != length(fields)
                    @printf(io, "| ")
                end
            end
            @printf(io, "\n")

            # print info for each frame
            for frame in reduced
                for (i, (l,f)) in enumerate(zip(field_lengths,fields))
                    @printf(io, "%-*s ", l, frame[f])
                    if i != length(fields)
                        @printf(io, "| ")
                    end
                end
                @printf(io, "\n")
            end

        end
            
        @showprogress for frame in reduced
            save(joinpath(plotting_folder, "frames_$(lpad(string(frame["FRAMENO"]), 4, '0')).png"), imview(frame, cmap=:matter))
        end

    end

end