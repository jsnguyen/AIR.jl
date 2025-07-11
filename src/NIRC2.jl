const _mask_dir = joinpath(@__DIR__, "..", "masks")

const NIRC2_bad_pixel_mask = begin
    maskfile = joinpath(_mask_dir, "bad_pixel_mask_20230101.fits")
    img = load(maskfile)            # from AstroImages or FITSIO
    Bool.(img.data)                 # convert to Bool matrix
end

# photoelectrons per adu
function get_NIRC2_gain(date_obs)
    d = Date(date_obs, dateformat"y-m-d")
    gain_date = Date(2023, 11, 20)
    if d < gain_date
        return 4.0  # pre-2023 gain
    else
        return 8.0  # post-2023 gain
    end
end

# readnoise in photoelectrons
function get_NIRC2_readnoise(sampmode)
    if sampmode == 2
        return 50.0
    elseif sampmode == 3
        return 15.0
    else
        return 0.0
    end
end