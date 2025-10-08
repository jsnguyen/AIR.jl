using Base.Threads
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Optim
using ADI
using Plots
using Serialization

using AIR
import AIR.crop 

function imshow(img; kwargs...)
    xs, ys = axes(img)
    heatmap(xs, ys, transpose(img); aspect_ratio=1,
            xlim=extrema(xs), ylim=extrema(ys), kwargs...)
end;

function adi()
    date = "2025-10-07"

    sequence_obslog_folder = "polmode/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "$(date)_split_sequences.toml")
    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = Obslog(sequence_obslog_path)

    key = "BD45598_dn_1"
    dn_frames = load(joinpath(sequence_obslog.paths.sequences_folder, "$(key)_aligned_frames.fits"), :)

    key = "BD45598_up_1"
    up_frames = load(joinpath(sequence_obslog.paths.sequences_folder, "$(key)_aligned_frames.fits"), :)

    combined = vcat(dn_frames, up_frames)

    sums = dn_frames .+ up_frames
    diffs = dn_frames .- up_frames

    angs = open(joinpath(sequence_obslog.paths.sequences_folder, "angles.jls"), "r") do io
        return deserialize(io)
    end

    angles = Float64[]
    for a in angs
        @info "ANGLE" a
        if a < 0
            a += 360
        end
        push!(angles, a)
    end


    cube = framelist_to_cube(sums)
    n_pca = [4,8,12,16,20]
    for np in n_pca
        alg = PCA(np)
        reduced = alg(cube, angles)
        save("reduced_sums_$(np).fits", reduced)
        #imshow(reduced)
        #savefig("reduced.png")
    end

    cube = framelist_to_cube(diffs)
    n_pca = [4,8,12,16,20]
    for np in n_pca
        alg = PCA(np)
        reduced = alg(cube, angles)
        save("reduced_diffs_$(np).fits", reduced)
        #imshow(reduced)
        #savefig("reduced.png")
    end

    cube = framelist_to_cube(dn_frames)
    n_pca = [4,8,12,16,20]
    for np in n_pca
        alg = PCA(np)
        reduced = alg(cube, angles)
        save("reduced_dn_$(np).fits", reduced)
        #imshow(reduced)
        #savefig("reduced.png")
    end

    cube = framelist_to_cube(up_frames)
    n_pca = [4,8,12,16,20]
    for np in n_pca
        alg = PCA(np)
        reduced = alg(cube, angles)
        save("reduced_up_$(np).fits", reduced)
        #imshow(reduced)
        #savefig("reduced.png")
    end



    alg = LOCI()
    reduced = alg(cube, angles)
    save("reduced_loci.fits", reduced)
    #imshow(reduced)
    #savefig("reduced.png")

    cube = framelist_to_cube(sums)
    #alg = alg(cube, angles)
    alg = PCA(4)
    psf = load(joinpath(sequence_obslog.paths.sequences_folder, "HD215806_dn_template_psf.fits"))
    fwhm = 5
    cc = contrast_curve(alg, cube, angles, psf; fwhm=fwhm)

    plot(
    cc.distance,
    [cc.contrast_corr cc.contrast],
    yscale=:log10,
    xlim=(0, NaN),
    label=["Student-t" "Gaussian"],
    ylabel="5-sigma contrast",
    xlabel="radius [px]"
    )

end

@autolog begin

    adi()

end