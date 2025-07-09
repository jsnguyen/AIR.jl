using OrderedCollections
using Glob
using FITSIO
using SkyCoords
using AstroAngles
using Statistics
using AstroImages

using AIR

function contains_target(data_folder, subfolder, coords; fk4_coords=nothing, n_cutoff=3, arcsec_threshold=100)

    filenames = glob("*.fits", joinpath(data_folder, subfolder))

    ra, dec = coords # assumed to be ICRS


    # if no FK4 coords, just use ICRS coordinates
    if fk4_coords !== nothing
        fk4_ra, fk4_dec = fk4_coords
    else
        fk4_ra, fk4_dec = ra, dec
    end



    target = String[]

    for fn in filenames
        frame = load(fn)

        conditions = [frame["CAMNAME"] == "narrow",
                      frame["GRSNAME"] == "clear",
                      !(haskey(frame, "SLITNAME") && (occursin("vortex", frame["SLITNAME"]) || occursin("corona", frame["SLITNAME"])))]

        if !all(conditions)
            @warn "Skipping file $(fn) due to conditions: $(conditions)"
            continue
        end

        if !(frame["RA"] isa Number) || !(frame["DEC"] isa Number)
            @warn "Skipping file $(fn) due to missing RA/DEC"
            continue
        end

        if frame["RADECSYS"] == "FK4"

            if fk4_coords === nothing
                @warn "Files are in FK4 but no Fk4 coordinates provided!"
            end

            @info "Using FK4 coordinates for $(fn)"
            ra = fk4_ra
            dec = fk4_dec
        else
            ra = ra
            dec = dec
        end

        radec_distance = deg2arcsec(small_angle_distance((ra, dec), (frame["RA"], frame["DEC"])))
        push!(distances, radec_distance)

        if radec_distance < arcsec_threshold && lowercase(frame["SHRNAME"]) == "open"

            @info "RADEC FRAME" filename=fn object=frame["OBJECT"] targname=frame["TARGNAME"] radec_distance=radec_distance
            push!(target, basename(fn))
        end
    end

    if length(target) <= n_cutoff
        println("Not enough science frames found, skipping obslog generation for $(data_folder)")
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
        println("No lamp-off flats found, assuming regular flats...")
        flats = flats_lampon
        flats_lampon = []
    end

    obslog = OrderedDict{String, Any}("data_folder" => data_folder,
                                      "subfolder" => subfolder,
                                      "date" => date,
                                      "sci" => sci,
                                      "flats" => flats,
                                      "flats_sky" => flats_sky,
                                      "flats_lampon" => flats_lampon,
                                      "flats_lampoff" => flats_lampoff,
                                      "darks" => darks)

    println("Found $(length(sci)) science frames")
    println("Found $(length(flats)) flats")
    println("Found $(length(flats_lampon)) lamp-on flats")
    println("Found $(length(flats_lampoff)) lamp-off flats")
    println("Found $(length(darks)) dark frames")

    toml_str = pretty_print_toml(obslog)

    open(obslog_filepath, "w") do io
        write(io, toml_str)
    end

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
        output_folder = "/Users/jsn/landing/projects/AIR.jl/reductions/obslogs"
        data_folder = joinpath(observations_folder, date)
        if contains_target(data_folder, subfolder, (ra, dec); fk4_coords=(fk4_ra, fk4_dec))
            if isdir(data_folder)
                generate_obs_logs(date, data_folder, subfolder, output_folder)
            end
        end
    end

end
