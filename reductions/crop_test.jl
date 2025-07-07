using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Plots

using AIR

function crop(img::AbstractArray, crop_size::Tuple{Int,Int}; center=nothing)
    h, w = size(img)
    crop_h, crop_w = crop_size

    if h==crop_h && w==crop_w
        return img  # No cropping needed, return the original image
    end

    if center === nothing
        center_row = Int(round(h/2))
        center_col = Int(round(w/2))
    else
        center_row, center_col = center
    end

    start_row = center_row - div(crop_h, 2)
    end_row   = start_row + crop_h - 1
    start_col = center_col - div(crop_w, 2)
    end_col   = start_col + crop_w - 1

    if start_row < 1 || start_col < 1 || end_row > h || end_col > w
        error("Crop would go out of bounds: img size $(size(img)), $(center_row), $(center_col), crop_size $(crop_size)")
    end

    img[start_row:end_row, start_col:end_col]
end

ones = ones((4,4))

heatmap(crop(ones, (2,2)))