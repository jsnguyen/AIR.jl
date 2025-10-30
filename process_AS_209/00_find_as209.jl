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

function is_lampoff_frame(frame, arcsec_threshold, lampoff_threshold)
    flat_el = 45.0 # degrees, elevation of the dome flat field frames, this seems to be consistent
    altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, frame["EL"])))

    lamp_flat_conditions = Dict{String,Any}("distance"       => altaz_distance <  arcsec_threshold,
                                            "dm loop status" => (lowercase(frame["WCDMSTAT"]) == "open" || lowercase(frame["WCDMSTAT"]) == "idle"),
                                            "dt loop status" => (lowercase(frame["WCDTSTAT"]) == "open" || lowercase(frame["WCDTSTAT"]) == "idle"),
                                            "shutter open"   => lowercase(frame["SHRNAME"]) == "open")

    if all(values(lamp_flat_conditions)) && median(frame) <= lampoff_threshold
        @info frame["FILENAME"] "is lamp off flat"
        return true
    end

    return false
end

function is_lampon_frame(frame, arcsec_threshold, lampoff_threshold)
    flat_el = 45.0 # degrees, elevation of the dome flat field frames, this seems to be consistent
    altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, frame["EL"])))

    lamp_flat_conditions = Dict{String,Any}("distance"       => altaz_distance <  arcsec_threshold,
                                            "dm loop status" => (lowercase(frame["WCDMSTAT"]) == "open" || lowercase(frame["WCDMSTAT"]) == "idle"),
                                            "dt loop status" => (lowercase(frame["WCDTSTAT"]) == "open" || lowercase(frame["WCDTSTAT"]) == "idle"),
                                            "shutter open"   => lowercase(frame["SHRNAME"]) == "open")


    if all(values(lamp_flat_conditions)) && median(frame) > lampoff_threshold
        @info frame["FILENAME"] "is lamp on flat"
        return true
    end

    return false
end

function is_sky_twilight_frame(frame)
    sky_twilight_flat_conditions = Dict{String,Any}("sky/twilight" => (occursin("sky", lowercase(frame["OBJECT"])) || occursin("twi", lowercase(frame["OBJECT"]))))

    if all(values(sky_twilight_flat_conditions))
        @info frame["FILENAME"] "is sky or twilight flat"
        return true
    end

    return false
end

function is_dark_frame(frame)
    dark_conditions = Dict{String,Any}("shutter closed" => lowercase(frame["SHRNAME"]) == "closed")

    if all(values(dark_conditions))
        @info frame["FILENAME"] "is dark frame"
        return true
    end

    return false
end

function sieve!(arr, original_indices, condition)
    indices = findall(condition, arr)
    actual_indices = original_indices[indices]
    deleteat!(original_indices, indices)
    deleteat!(arr, indices)
    return actual_indices
end

function sort_frames(frames, filenames; lampoff_threshold=100.0, arcsec_threshold=100.0)

    required_header_keywords = ["EL", "WCDMSTAT", "WCDTSTAT", "OBJECT", "SHRNAME"]

    for (i, f) in enumerate(frames)
        for key in required_header_keywords
            if !haskey(f, key)
                @warn "Frame $(filenames[i]) is missing required header keyword: $key, removing from list..."
                deleteat!(frames, i)
                deleteat!(filenames, i)
            end
        end
    end

    # start with all the indices then "sieve" out the different types

    original_indices = collect(eachindex(frames))

    flats_lampon_indices = sieve!(frames, original_indices, x -> is_lampon_frame(x, arcsec_threshold, lampoff_threshold))
    flats_lampoff_indices = sieve!(frames, original_indices, x -> is_lampoff_frame(x, arcsec_threshold, lampoff_threshold))
    flats_sky_indices = sieve!(frames, original_indices, x -> is_sky_twilight_frame(x))
    darks_indices = sieve!(frames, original_indices, x -> is_dark_frame(x))
    sci_indices = original_indices

    flats_lampon = filenames[flats_lampon_indices]
    flats_lampoff = filenames[flats_lampoff_indices]
    flats_sky = filenames[flats_sky_indices]
    darks = filenames[darks_indices]
    sci = filenames[sci_indices]

    if isempty(flats_lampoff)
        @info "No lamp-off flats found, assuming regular flats..."
        flats = flats_lampon
        flats_lampon = String[]
    else
        flats = String[]
    end

    @info "Found $(length(sci)) science frames"
    @info "Found $(length(flats)) flats"
    @info "Found $(length(flats_lampon)) lamp-on flats"
    @info "Found $(length(flats_lampoff)) lamp-off flats"
    @info "Found $(length(darks)) dark frames"

    return sci, flats, flats_sky, flats_lampon, flats_lampoff, darks
end

@stage function find_as209()

    ra      = hms"16 49 15.3034917000"deg
    dec     = dms"-14 22 08.643317664"deg

    fk4_ra  = hms"16 46 25.3250542612"deg
    fk4_dec = dms"-14 16 57.045137674"deg

    @info "=== AS 209 coordinates ==="
    @info "ICRS RA DEC" ra dec
    @info "FK4 RA DEC" fk4_ra fk4_dec

    dates_containing_as209 = String[]
    for date in readdir(context["observations_folder"])

        paths = ObslogPaths(date, joinpath(context["observations_folder"], date))
        @info "Processing date" date=date

        frames = AstroImage[]
        filenames = String[]
        for fn in glob("*.fits", joinpath(paths.data_folder, "raw"))
            push!(frames, load(fn))
            push!(filenames, fn)
        end

        make_frametable(frames, joinpath(paths.data_folder, "$(date)_frames_table.txt"))

        if contains_target(frames, filenames, (ra, dec); fk4_coords=(fk4_ra, fk4_dec))
            push!(dates_containing_as209, date)
        end

    end

    @info "Dates containing AS 209:" dates=dates_containing_as209
    
    return dates_containing_as209

end