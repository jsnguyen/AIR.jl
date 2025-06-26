using OrderedCollections

using Glob
using FITSIO
using SkyCoords
using AstroAngles
using Plots
using Statistics

using AIR


function generate_obs_logs(date, data_folder, subfolder, output_folder, lampoff_threshold=100.0)

    if !isdir(output_folder)
        mkpath(output_folder)
    end

    obslog_filepath = joinpath(output_folder, "$(date)_AS_209.toml")

    filenames = glob("*.fits", joinpath(data_folder, subfolder))

    icrs_ra = hms"16 49 15.3034917000"deg
    icrs_dec = dms"-14 22 08.643317664"deg

    fk4_ra = hms"16 46 25.3250542612"deg
    fk4_dec = dms"-14 16 57.045137674"deg

    arcsec_threshold = 100 # arcsec
    @info "=== AS 209 coordinates ==="
    @info "ICRS RA DEC" icrs_ra icrs_dec
    @info "FK4 RA DEC" fk4_ra fk4_dec

    flat_az = 170.0 # this doesnt seem needed? az not necessarily at 170 deg in early Keck data
    flat_el = 45.0

    small_angle_distance = ((ra_a, dec_a), (ra_b, dec_b)) -> sqrt(((ra_a-ra_b)*cos(deg2rad(dec_a)))^2 + (dec_a-dec_b)^2) # units of degrees

    sci = []
    flats = []
    flats_sky = []
    flats_lampon = []
    flats_lampoff = []
    darks = []

    frames = []
    headers = []

    distances = []
    for fn in filenames
        hdu = FITS(fn, "r")
        d = read(hdu[1])
        h = read_header(hdu[1])
        push!(frames, d)
        push!(headers, h)

        conditions = [h["CAMNAME"] == "narrow",
                      h["GRSNAME"] == "clear",
                      !(haskey(h, "SLITNAME") && (occursin("vortex", h["SLITNAME"]) || occursin("corona", h["SLITNAME"])))]

        if !all(conditions)
            @warn "Skipping file $(fn) due to conditions: $(conditions)"
            continue
        end

        if !(h["RA"] isa Number) || !(h["DEC"] isa Number)
            @warn "Skipping file $(fn) due to missing RA/DEC"
            continue
        end

        if h["RADECSYS"] == "FK4"
            @info "Using FK4 coordinates for $(fn)"
            ra = fk4_ra
            dec = fk4_dec
        else
            ra = icrs_ra
            dec = icrs_dec
        end

        radec_distance = deg2arcsec(small_angle_distance((ra, dec), (h["RA"], h["DEC"])))
        altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, h["EL"])))
        push!(distances, radec_distance)

        if radec_distance < arcsec_threshold && lowercase(h["SHRNAME"]) == "open"

            @info "RADEC FRAME" filename=fn object=h["OBJECT"] targname=h["TARGNAME"] radec_distance=radec_distance
            push!(sci, basename(fn))

        elseif altaz_distance <  arcsec_threshold && lowercase(h["WCDMSTAT"]) == "open" && lowercase(h["WCDTSTAT"]) == "open"

            if median(d) < lampoff_threshold
                @info "LAMPOFF FRAME" filename=fn object=h["OBJECT"] targname=h["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampoff, basename(fn))
            else
                @info "LAMPON FRAME" filename=fn object=h["OBJECT"] targname=h["TARGNAME"] altaz_distance=altaz_distance
                push!(flats_lampon, basename(fn))
            end

        elseif occursin("sky", lowercase(h["OBJECT"]))

            @info "SKY FLAT FRAME" filename=fn object=h["OBJECT"] targname=h["TARGNAME"]
            push!(flats_sky, basename(fn))


        elseif h["SHRNAME"] == "closed"

            @info "DARK FRAME" filename=fn object=h["OBJECT"] targname=h["TARGNAME"]
            push!(darks, basename(fn))

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

    if length(sci) <= 3
        println("Not enough science frames found, skipping obslog generation for $(date)")
        return
    end

    toml_str = pretty_print_toml(obslog)

    open(obslog_filepath, "w") do io
        write(io, toml_str)
    end

    return distances

end

autolog("$(@__FILE__).log") do

    observations_folder = "/Users/jsn/landing/projects/AIR.jl/data/"

    subfolder = "raw"
    for date in readdir(observations_folder)
        output_folder = "/Users/jsn/landing/projects/AIR.jl/reductions/obslogs"
        data_folder = joinpath(observations_folder, date)
        if isdir(data_folder)
            distances = generate_obs_logs(date, data_folder, subfolder, output_folder)
        end
    end
end
