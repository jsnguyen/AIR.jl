using Printf
using OrderedCollections
using Glob
using FITSIO
using SkyCoords
using AstroAngles
using Statistics
using AstroImages

using AIR

function make_table(date, data_folder, subfolder)

    filenames = glob("*.fits", joinpath(data_folder, subfolder))

    fields = ["FILENAME", "OBJECT", "TARGNAME", "RA", "DEC", "CAMNAME", "DATE-OBS", "UTC", "ITIME", "COADDS", "FILTER", "CAMNAME", "GRSNAME", "SLITNAME", "NAXIS1", "NAXIS2"]
    field_lengths = [18, 16, 16, 12, 12, 8, 10, 11, 8, 8, 20, 8, 8, 10, 6, 6]

    table_filename = joinpath(data_folder, "$(date)_frames_table.txt")
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
        for fn in filenames
            frame = load(fn)
            for (i, (l,f)) in enumerate(zip(field_lengths,fields))
                if !haskey(frame, f)
                    @printf(io, "%-*s ", l, "")
                    continue
                else
                    @printf(io, "%-*s ", l, frame[f])
                end

                if i != length(fields)
                    @printf(io, "| ")
                end

            end
            @printf(io, "\n")
        end

    end
end

function contains_target(data_folder, subfolder, coords; fk4_coords=nothing, n_cutoff=3, arcsec_threshold=100)

    filenames = glob("*.fits", joinpath(data_folder, subfolder))

    target = String[]

    for fn in filenames
        frame = load(fn)

        conditions = Dict{String,Any}("NARROWCAM" => frame["CAMNAME"] == "narrow",
                                      "CLEAR GRS" => frame["GRSNAME"] == "clear")
                                      #"OPEN SLIT" => !(haskey(frame, "SLITNAME") && (occursin("vortex", frame["SLITNAME"]) || occursin("corona", frame["SLITNAME"]))))

        if !all(values(conditions))
            @warn "Skipping file $(fn) due to conditions" conditions=conditions
            continue
        end

        if !(frame["RA"] isa Number) || !(frame["DEC"] isa Number)
            @warn "Skipping file $(fn) due to missing RA/DEC"
            continue
        end

        # need this here or they will get overwritten on subsequent loops
        # for some reason, the coordinate system can vary within an epoch
        ra, dec = coords

        if frame["RADECSYS"] == "FK4"
            @info "Using FK4 coordinates for $(fn)"
            if fk4_coords !== nothing
                ra, dec = fk4_coords
            else
                @warn "Files are in FK4 but no Fk4 coordinates provided!"
            end
        end

        radec_distance = deg2arcsec(small_angle_distance((ra, dec), (frame["RA"], frame["DEC"])))
        @info "Coordinates" object=frame["OBJECT"] target=frame["TARGNAME"] ra=ra dec=dec frame_ra=frame["RA"] frame_dec=frame["DEC"] radec_distance=radec_distance radecsys=frame["RADECSYS"]

        if radec_distance < arcsec_threshold && lowercase(frame["SHRNAME"]) == "open"

            @info "RADEC FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"] radec_distance=radec_distance
            push!(target, basename(fn))
        end
    end

    if length(target) <= n_cutoff
        @warn "Not enough science frames found (only found $(length(target))), skipping obslog generation for $(data_folder)"
        return false
    end

    return true

end

function generate_obs_logs(date, data_folder, subfolder, output_folder; lampoff_threshold=100.0, arcsec_threshold=100.0)

    if !isdir(output_folder)
        mkpath(output_folder)
    end

    obslog_filepath = joinpath(output_folder, "$(date)_obslog.toml")
    filenames = glob("*.fits", joinpath(data_folder, subfolder))

    flat_el = 45.0 # degrees, elevation of the flat field frames

    sci = []
    flats = []
    flats_sky = []
    flats_lampon = []
    flats_lampoff = []
    darks = []

    for fn in filenames
        frame = load(fn)

        altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, frame["EL"])))

        if altaz_distance <  arcsec_threshold && lowercase(frame["WCDMSTAT"]) == "open" && lowercase(frame["WCDTSTAT"]) == "open"

            if median(frame) < lampoff_threshold
                @info "LAMPOFF FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampoff, basename(fn))
            else
                @info "LAMPON FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampon, basename(fn))
            end

        elseif occursin("sky", lowercase(frame["OBJECT"])) | occursin("twi", lowercase(frame["OBJECT"]))

            @info "SKY FLAT FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(flats_sky, basename(fn))

        elseif frame["SHRNAME"] == "closed"

            @info "DARK FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(darks, basename(fn))

        else

            @info "SCI FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(sci, basename(fn))

        end

    end

    if flats_lampoff == []
        @info "No lamp-off flats found, assuming regular flats..."
        flats = flats_lampon
        flats_lampon = []
    end

    obslog = OrderedDict{String, Any}("data_folder" => data_folder,
                                      "date" => date,
                                      subfolder => OrderedDict{String, Any}(
                                        "sci" => sci,
                                        "flats" => flats,
                                        "flats_sky" => flats_sky,
                                        "flats_lampon" => flats_lampon,
                                        "flats_lampoff" => flats_lampoff,
                                        "darks" => darks
                                      ))


    @info "Found $(length(sci)) science frames"
    @info "Found $(length(flats)) flats"
    @info "Found $(length(flats_lampon)) lamp-on flats"
    @info "Found $(length(flats_lampoff)) lamp-off flats"
    @info "Found $(length(darks)) dark frames"

    toml_str = pretty_print_toml(obslog)

    open(obslog_filepath, "w") do io
        write(io, toml_str)
    end

    return obslog

end

@autolog begin

    observations_folder = "/Users/jsn/landing/projects/AIR.jl/data/"

    ra      = hms"16 49 15.3034917000"deg
    dec     = dms"-14 22 08.643317664"deg

    fk4_ra  = hms"16 46 25.3250542612"deg
    fk4_dec = dms"-14 16 57.045137674"deg

    @info "=== AS 209 coordinates ==="
    @info "ICRS RA DEC" ra dec
    @info "FK4 RA DEC" fk4_ra fk4_dec

    subfolder = "raw"

    for date in readdir(observations_folder)
        if !isdir(joinpath(observations_folder, date))
            continue
        end
        output_folder = "/Users/jsn/landing/projects/AIR.jl/pipeline/obslogs"
        data_folder = joinpath(observations_folder, date)

        make_table(date, data_folder, subfolder)

        if contains_target(data_folder, subfolder, (ra, dec); fk4_coords=(fk4_ra, fk4_dec))
            if isdir(data_folder)
                generate_obs_logs(date, data_folder, subfolder, output_folder)
            end
        end

    end

end
