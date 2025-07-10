using Glob
using Printf
using ProgressMeter
using AstroImages

using AIR

function make_table(obslog, obslog_folder)
    fields = ["RED-FN", "OBJECT", "TARGNAME", "RA", "DEC", "CAMNAME", "DATE-OBS", "UTC", "ITIME", "COADDS", "FILTER", "NAXIS1", "NAXIS2"]
    field_lengths = [18, 16, 16, 12, 12, 8, 10, 11, 8, 8, 20, 6, 6]

    table_filename = joinpath(obslog_folder, "$(obslog["date"])_reduced_frames_table.txt")
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
        for frame in obslog.reduced_sci
            for (i, (l,f)) in enumerate(zip(field_lengths,fields))
                @printf(io, "%-*s ", l, frame[f])
                if i != length(fields)
                    @printf(io, "| ")
                end
            end
            @printf(io, "\n")
        end

    end
end

@autolog begin

    obslog_folder = "pipeline/obslogs"
    for obslog_filename in Glob.glob(joinpath(obslog_folder,"*_reduced.toml"))

        @info "Loading obslog from" obslog_filename

        obslog = Obslog(obslog_filename)
        reduced = obslog.reduced_sci

        #
        # Plot
        #

        plotting_folder = joinpath(obslog["data_folder"], "plots")
        make_and_clear(plotting_folder, "frames_*.png")

        make_table(obslog, obslog_folder)
            
        Threads.@threads for i in eachindex(reduced)
            frame = reduced[i]
            save_filename = joinpath(plotting_folder, "frames_$(lpad(string(frame["FRAMENO"]), 4, '0')).png")
            @info "Saving frame to $(save_filename)"
            save(save_filename, imview(frame, cmap=:matter))
        end

    end

end