using Glob
using Printf
using ProgressMeter
using AstroImages
using Distributed

using AIR

@autolog begin

    obslog_folder = "pipeline/obslogs"
    for obslog_filename in Glob.glob(joinpath(obslog_folder,"*_reduced.toml"))

        @info "Loading obslog from" obslog_filename
        obslog = Obslog(obslog_filename)

        fields = ["FILENAME", "OBJECT", "TARGNAME", "RA", "DEC", "CAMNAME", "DATE-OBS", "UTC", "ITIME", "COADDS", "FILTER", "NAXIS1", "NAXIS2"]
        make_frametable(obslog.reduced_sci, obslog.paths.table_file; fields=fields)

        #
        # Plot
        #
        
        make_and_clear(obslog.paths.plots_folder, "frames_*.png")

        p = Progress(length(obslog.reduced_sci))
        counter = Threads.Atomic{Int}(0)

        Threads.@threads for i in eachindex(obslog.reduced_sci)
            frame = obslog.reduced_sci[i]
            save_filename = joinpath(obslog.paths.plots_folder, "frames_$(lpad(string(frame["FRAMENO"]), 4, '0')).png")
            save(save_filename, imview(frame, cmap=:matter))

            # no thread locking here for SPEED
            Threads.atomic_add!(counter, 1)
            update!(p, counter[])
        end

    end

end