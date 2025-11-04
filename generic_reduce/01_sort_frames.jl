using Printf
using OrderedCollections
using Glob
using FITSIO
using SkyCoords
using AstroAngles
using Statistics
using AstroImages

using AIR

function is_lampoff_frame(frame, arcsec_threshold, lampoff_threshold)
    flat_el = 45.0 # degrees, elevation of the dome flat field frames, this seems to be consistent
    altaz_distance = deg2arcsec(small_angle_distance((0, flat_el), (0, frame["EL"])))

    lamp_flat_conditions = Dict{String,Any}("distance"       => altaz_distance <  arcsec_threshold,
                                            "dm loop status" => (lowercase(frame["WCDMSTAT"]) == "open" || lowercase(frame["WCDMSTAT"]) == "idle"),
                                            "dt loop status" => (lowercase(frame["WCDTSTAT"]) == "open" || lowercase(frame["WCDTSTAT"]) == "idle"),
                                            "shutter open"   => lowercase(frame["SHRNAME"]) == "open")

    if all(values(lamp_flat_conditions)) && median(frame) <= lampoff_threshold
        @info "$(frame["FILENAME"]) is lamp off flat"
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
        @info "$(frame["FILENAME"]) is lamp on flat"
        return true
    end

    return false
end

function is_sky_twilight_frame(frame)
    sky_twilight_flat_conditions = Dict{String,Any}("sky/twilight" => (occursin("sky", lowercase(frame["OBJECT"])) || occursin("twi", lowercase(frame["OBJECT"]))))

    if all(values(sky_twilight_flat_conditions))
        @info "$(frame["FILENAME"]) is sky or twilight flat"
        return true
    end

    return false
end

function is_dark_frame(frame)
    dark_conditions = Dict{String,Any}("shutter closed" => lowercase(frame["SHRNAME"]) == "closed")

    if all(values(dark_conditions))
        @info "$(frame["FILENAME"]) is dark frame"
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

@stage function sort_frames(; lampoff_threshold=100.0, arcsec_threshold=100.0)

    data_folder = context["paths"].data_folder
    date = context["paths"].date

    frames = AstroImage[]
    filenames = String[]
    for fn in glob("*.fits", joinpath(data_folder, "raw"))
        push!(frames, load(fn))
        push!(filenames, fn)
    end

    make_frametable(frames, joinpath(data_folder, "$(date)_frames_table.txt"))

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