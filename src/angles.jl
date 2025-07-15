"""
Calculate and return the starting angle, angular smear, and mean angle to north
of an image from NIRC2 narrow cam, given its headers.

The calculation is adapted from pyKLIP.
"""

float_or_parse_hexages(num::Number) = num
float_or_parse_hexages(str::AbstractString) = AstroLib.ten(str)

function calculate_north_angle(headers)
    zp_offset = -0.262 # From Service et al 2016
    rotator_mode = headers["ROTMODE"]
    rotator_position = headers["ROTPOSN"] # Degrees
    instrument_angle = headers["INSTANGL"] # Degrees

    pa_deg = 0.0
    if rotator_mode == "vertical angle"
        if haskey(headers, "PARANTEL")
            parang = headers["PARANTEL"]
        else
            parang = headers["PARANG"]
        end
        pa_deg = parang + rotator_position - instrument_angle + zp_offset

    # fixed bug here where position angle mode was not being handled correctly
    # don't need all the fancy smear stuff if we're in position angle mode
    elseif rotator_mode == "position angle"
        pa_deg = rotator_position - instrument_angle + zp_offset
        return pa_deg, [pa_deg]
    elseif rotator_mode == "stationary"
        # TODO: Handle case where the instrument rotator is stationary.
        return NaN, [NaN]
    else
        throw(ArgumentError("Unknown rotator mode " * rotator_mode))
    end

    # Keck rotator bug.
            # parang = headers["PARANG"] # Degrees
    # This does not appear to actually work
    diff = 0.0
    # if haskey(headers, "ROTNORTH")
    #     pa_deg_idl = -headers["ROTNORTH"]
    #     # diff += (pa_deg_idl - pa_deg)
    #     # println(diff)
    #     @info "Using ROTNORTH header from IDL pipeline" maxlog=1
    # end

    # Now calculate smear

    # Get info for PA smearing calculation.
    # epochobj = headers["DATE-OBS"]
    # name = headers["TARGNAME"]
    expref = headers["ITIME"]
    coaddref = headers["COADDS"]
    sampref = headers["SAMPMODE"]
    msrref = headers["MULTISAM"]
    xdimref = headers["NAXIS1"]
    # ydimref = headers["NAXIS2"]
    dec = float_or_parse_hexages(headers["DEC"]) + headers["DECOFF"]

    if haskey(headers, "TOTEXP")
        totexp = headers["TOTEXP"]
    else
        # Calculate total time of exposure (integration + readout).
        if sampref == 2
            totexp = (expref + 0.18 * (xdimref / 1024.0)^2) * coaddref
        end
        if sampref == 3
            totexp =
                (expref + (msrref - 1) * 0.18 * (xdimref / 1024.0)^2) * coaddref
        end
    end
    # tinteg = totexp # [seconds]
    totexp = totexp / 3600.0 # [hours]

    # Get hour angle at start of exposure.
    tmpahinit = float_or_parse_hexages(headers["HA"]) # [deg]
    ahobs = 24.0 * tmpahinit / 360.0 # [hours]

    if totexp * 3600 > 1 # If greater than 1 second...
        # Estimate vertical position angle at each second of the exposure.
        vp = Float64[]
        vpref = 0.0
        for j in 0:(3600*totexp-1)
            ahtmp = ahobs + (j + 1.0 + 0.001) / 3600.0 # hours
            # TODO: par_angle and observatory_latitude
            push!(vp, par_angle(ahtmp, dec, observatory_lat))
            if j == 0
                vpref = vp[1]
            end
        end

        # Handle case where PA crosses 0 <--> 360.
        vp[vp.<0] .+= 360
        vp[vp.>360.0] .+= 360

        if vpref < 0
            vpref += 360
        end
        if vpref > 360
            vpref -= 360
        end

        # TODO: This is some crazy code... Should use use meandegrees()
        # that uses atan() to avoid this.
        # Check that images near PA=0 are handled correctly.
        if any(vp .> 350) && any(vp .< 10)
            vp[vp.>350] .-= 360
        end
        vpmean = mean(angle for angle = vp if isfinite(angle))

        if vpmean < 60 && vpref > 350
            vpmean += 360
            vp .+= 360
        end
        pa_deg_mean = pa_deg + (vpmean - vpref)

        # angle_start = pa_deg
        # angle_smear = vpmean - vpref
        angle_mean = pa_deg_mean
        return angle_mean+diff, vp.+diff
    else # Total exposure less than one second - treat as no smear
        return pa_deg+diff, [pa_deg+diff]
    end
end

"""
Compute the parallactic angle, given hour angle (HA [hours]),
declination (dec [deg]), and latitude (lat [deg]).  Returns
parallactic angle in [deg].
Source: pyKLIP
"""
function par_angle(HA, dec, lat)

    HA_rad = deg2rad(HA * 15.0) # [hours] -> [rad]
    dec_rad = deg2rad(dec)   # [deg] -> [rad]
    lat_rad = deg2rad(lat)   # [deg] -> [rad]

    parallang =
        -atan(
            -sin(HA_rad),  # [rad]
            cos(dec_rad) * tan(lat_rad) - sin(dec_rad) * cos(HA_rad),
        )

    return rad2deg(parallang) # [deg]
end

"""
    rotate_image_center(img::AstroImage, angle_degrees; fillval=0.0)

Rotate an image about its center by the specified angle in degrees.
"""
function rotate_image_center(img::AbstractArray, angle_degrees; fillval=0.0)
    angle_rad = deg2rad(angle_degrees)
    
    rows, cols = size(img)
    center = rows/2 + 0.5, cols/2 + 0.5
    
    # Create rotation about specified center point
    rotation = recenter(RotMatrix(angle_rad), center)
    
    rotated_data = warp(img, rotation, axes(img), fill=fillval)
    
    return rotated_data
end

function rotate_image_center(img::AstroImage, angle_degrees; fillval=0.0)
    angle_rad = deg2rad(angle_degrees)

    
    rows, cols = size(img)
    center = rows/2 + 0.5, cols/2 + 0.5
    
    # Create rotation about specified center point
    rotation = recenter(RotMatrix(angle_rad), center)
    
    rotated_data = warp(img.data, rotation, axes(img.data), fill=fillval)
    
    return AstroImage(rotated_data, img.header)
end