using Printf
using OrderedCollections
using Glob
using FITSIO
using SkyCoords
using AstroAngles
using Statistics
using AstroImages

using AIR


function contains_target(frames, filenames, coords; fk4_coords=nothing, n_cutoff=3, arcsec_threshold=100)

    n_target_frames = 0
    for (i, frame) in enumerate(frames)

        conditions = Dict{String,Any}("NARROWCAM" => frame["CAMNAME"] == "narrow",
                                      "CLEAR GRS" => frame["GRSNAME"] == "clear")
                                      #"OPEN SLIT" => !(haskey(frame, "SLITNAME") && (occursin("vortex", frame["SLITNAME"]) || occursin("corona", frame["SLITNAME"]))))

        if !all(values(conditions))
            @warn "Skipping file $(filenames[i]) due to conditions" conditions=conditions
            continue
        end

        if !(frame["RA"] isa Number) || !(frame["DEC"] isa Number)
            @warn "Skipping file $(filenames[i]) due to missing RA/DEC"
            continue
        end

        # need this here or they will get overwritten on subsequent loops
        # for some reason, the coordinate system can vary within an epoch
        ra, dec = coords

        if frame["RADECSYS"] == "FK4"
            @info "Using FK4 coordinates for $(filenames[i])"
            if fk4_coords !== nothing
                ra, dec = fk4_coords
            else
                @warn "Files are in FK4 but no Fk4 coordinates provided!"
            end
        end

        radec_distance = deg2arcsec(small_angle_distance((ra, dec), (frame["RA"], frame["DEC"])))
        @info "Coordinates" object=frame["OBJECT"] target=frame["TARGNAME"] ra=ra dec=dec frame_ra=frame["RA"] frame_dec=frame["DEC"] radec_distance=radec_distance radecsys=frame["RADECSYS"]

        if radec_distance < arcsec_threshold && lowercase(frame["SHRNAME"]) == "open"

            @info "RADEC FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"] radec_distance=radec_distance
            n_target_frames += 1
        end
    end

    if n_target_frames <= n_cutoff
        @warn "Not enough science frames found (only found $(n_target_frames)), skipping obslog generation"
        return false
    end

    return true

end

function sort_frames(frames, filenames; lampoff_threshold=100.0, arcsec_threshold=100.0)

    flat_el = 45.0 # degrees, elevation of the flat field frames

    sci = []
    flats = []
    flats_sky = []
    flats_lampon = []
    flats_lampoff = []
    darks = []

    for (i, frame) in enumerate(frames)

        altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, frame["EL"])))

        if altaz_distance <  arcsec_threshold && lowercase(frame["WCDMSTAT"]) == "open" && lowercase(frame["WCDTSTAT"]) == "open"

            if median(frame) < lampoff_threshold
                @info "LAMPOFF FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampoff, filenames[i])
            else
                @info "LAMPON FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampon, filenames[i])
            end

        elseif occursin("sky", lowercase(frame["OBJECT"])) | occursin("twi", lowercase(frame["OBJECT"]))

            @info "SKY FLAT FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(flats_sky, filenames[i])

        elseif frame["SHRNAME"] == "closed"

            @info "DARK FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(darks, filenames[i])

        else

            @info "SCI FRAME" filename=filenames[i] object=frame["OBJECT"] targname=frame["TARGNAME"]
            push!(sci, filenames[i])

        end

    end

    if flats_lampoff == []
        @info "No lamp-off flats found, assuming regular flats..."
        flats = flats_lampon
        flats_lampon = []
    end

    @info "Found $(length(sci)) science frames"
    @info "Found $(length(flats)) flats"
    @info "Found $(length(flats_lampon)) lamp-on flats"
    @info "Found $(length(flats_lampoff)) lamp-off flats"
    @info "Found $(length(darks)) dark frames"

    return sci, flats, flats_sky, flats_lampon, flats_lampoff, darks
end

function make_obslog(data_folder, date, obslog_filepath, sci, flats, flats_sky, flats_lampon, flats_lampoff, darks) 

    obslog = OrderedDict{String, Any}("data_folder" => data_folder,
                                      "date" => date,
                                      "raw" => OrderedDict{String, Any}(
                                        "sci" => sci,
                                        "flats" => flats,
                                        "flats_sky" => flats_sky,
                                        "flats_lampon" => flats_lampon,
                                        "flats_lampoff" => flats_lampoff,
                                        "darks" => darks
                                      ))

    write_toml(obslog_filepath, obslog)

end

function generate_obslogs_as209()
    observations_folder = "/Users/jsn/landing/projects/AIR.jl/data/"

    ra      = hms"16 49 15.3034917000"deg
    dec     = dms"-14 22 08.643317664"deg

    fk4_ra  = hms"16 46 25.3250542612"deg
    fk4_dec = dms"-14 16 57.045137674"deg

    @info "=== AS 209 coordinates ==="
    @info "ICRS RA DEC" ra dec
    @info "FK4 RA DEC" fk4_ra fk4_dec

    for date in readdir(observations_folder)

        if !isdir(joinpath(observations_folder, date))
            continue
        end

        output_folder = "/Users/jsn/landing/projects/AIR.jl/pipeline/obslogs"
        mkpath(output_folder)

        data_folder = joinpath(observations_folder, date)

        frames = AstroImage[]
        filenames = String[]
        for fn in glob("*.fits", joinpath(data_folder, "raw"))
            push!(frames, load(fn))
            push!(filenames, fn)
        end

        make_frametable(frames, joinpath(data_folder, "$(date)_frames_table.txt"))

        if contains_target(frames, filenames, (ra, dec); fk4_coords=(fk4_ra, fk4_dec))
            if isdir(data_folder)
                sorted_frames = sort_frames(frames, filenames)
                make_obslog(data_folder, date, joinpath(output_folder, "$(date)_obslog.toml"), sorted_frames...)
            end
        end

    end
end

@autolog begin

    generate_obslogs_as209()

end
