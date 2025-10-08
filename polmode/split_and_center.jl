using Printf

using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations

using AIR
import AIR.crop 

function split_image_save(img::AstroImage, out_dir::AbstractString, basename::AbstractString; direction::Symbol=:vertical, mid::Union{Nothing,Int}=nothing)
    """
    Split `img` in half and save both halves as FITS.

    - direction=:vertical (splits left/right across the vertical axis) or :horizontal (splits up/dn across the horizontal axis)
    - out_dir: directory to save files
    - basename: base filename (no extension)
    - mid: optional split index (column for vertical, row for horizontal)
    Returns tuple of saved file paths.
    """
    mkpath(out_dir)

    data = img.data
    nrows, ncols = size(data)

    if mid === nothing
        mid = direction === :vertical ? fld(ncols, 2) : fld(nrows, 2)
    end
    mid = Int(mid)

    if direction === :vertical
        # split across vertical axis -> left / right (l / r)
        left_data = copy(data[1:mid, :])
        right_data = copy(data[mid+1:nrows, :])
        left  = AstroImage(left_data, deepcopy(img.header))
        right = AstroImage(right_data, deepcopy(img.header))
        left["SPLIT_PART"]  = "l"
        right["SPLIT_PART"] = "r"
        left_path  = joinpath(out_dir, "$(basename)_l.fits")
        right_path = joinpath(out_dir, "$(basename)_r.fits")
        save(left_path, left)
        save(right_path, right)
        return left_path, right_path

    elseif direction === :horizontal
        # split across horizontal axis -> up / dn
        top_data = copy(data[:, mid+1:ncols])
        bottom_data  = copy(data[:, 1:mid])
        top    = AstroImage(top_data, deepcopy(img.header))
        bottom = AstroImage(bottom_data, deepcopy(img.header))
        top["SPLIT_PART"]    = "up"
        bottom["SPLIT_PART"] = "dn"
        top_path    = joinpath(out_dir, "$(basename)_up.fits")
        bottom_path = joinpath(out_dir, "$(basename)_dn.fits")
        save(top_path, top)
        save(bottom_path, bottom)
        return top_path, bottom_path

    else
        error("direction must be :vertical or :horizontal")
    end
end

function split_epochs()

    date = "2025-10-07"

    sequence_obslog_folder = "live_ingestion/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "$(date)_sequences.toml")
    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    split_folder = joinpath(sequence_obslog.paths.data_folder, "reduced")
    @info split_folder
    mkpath(split_folder)

    for key in keys(sequence_obslog.sequences)

        @info "Processing sequence" key

        frames = load_frames(sequence_obslog, key) 

        for f in frames
            split_image_save(f, split_folder, "$(chop(f["RED-FN"], tail=5))"; direction=:horizontal)
        end
    end

end

@autolog begin
    split_epochs()
end