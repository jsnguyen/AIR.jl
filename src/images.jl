
function make_sigma_clip_mask(image_data::AbstractMatrix, n_sigma::Real = 9.0)
    if isempty(image_data)
        return BitMatrix(zeros(Bool, 0, 0))
    end

    mean_val = median(image_data)
    std_dev = std(image_data)
    threshold = mean_val + (n_sigma * std_dev)
    mask = image_data .> threshold

    return mask
end
