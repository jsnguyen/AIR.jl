using Glob
using Printf
using ProgressMeter
using AstroImages
using Distributed

using AIR

@stage function make_table_and_plots(reduced_sci)
    paths = context["paths"]

    fields = ["FILENAME", "OBJECT", "TARGNAME", "RA", "DEC", "CAMNAME", "DATE-OBS", "UTC", "ITIME", "COADDS", "SLITNAME", "AZ", "EL", "PARANG", "FILTER", "NAXIS1", "NAXIS2"]
    make_frametable(reduced_sci, paths.table_file; fields=fields)

    make_and_clear(paths.plots_folder, "frames_*.png")

    p = Progress(length(reduced_sci))
    counter = Threads.Atomic{Int}(0)

    Threads.@threads for i in eachindex(reduced_sci)
        frame = reduced_sci[i]
        save_filename = joinpath(paths.plots_folder, "frames_$(lpad(string(frame["FRAMENO"]), 4, '0')).png")
        save(save_filename, imview(frame, cmap=:matter))

        # no thread locking here for SPEED
        Threads.atomic_add!(counter, 1)
        update!(p, counter[])
    end

    return

end