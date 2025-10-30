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

function do_pca(cube, angles, paths, n_pca)

    algs = ADI.ADIAlgorithm[]
    reduced = AstroImage[]
    for n in n_pca
        @info "Calculating PCA reduction for PCA n = $(n)..."

        alg = PCA(n)
        res = AstroImage(alg(cube, angles))

        res["N_PCA"] = n

        push!(algs, alg)
        push!(reduced, res)

        save(joinpath(paths.sequences_folder, "reduced_pca_$(n).fits"), res)

        xs, ys = axes(res)
        heatmap(xs, ys, transpose(res); aspect_ratio=1, xlim=extrema(xs), ylim=extrema(ys), clims=(-5,5))
        savefig(joinpath(paths.sequences_folder, "reduced_pca_$(n).png"))
    end

    return algs, reduced

end

function calc_contrast_curves(algs, cube, angles, template_psf, fwhm, paths, n_pca)

    for (i, n) in enumerate(n_pca)
        @info "Calculating contrast curve for PCA n = $(n)..."
        cc = contrast_curve(algs[i], cube, angles, template_psf; fwhm=fwhm)
        serialize(joinpath(paths.sequences_folder, "contrast_curve_pca_$(n).jls"), cc)
    end

end

function plot_contrast_curves(paths, n_pca)

    plot()

    for n in n_pca
        cc = deserialize(joinpath(paths.sequences_folder, "contrast_curve_pca_$(n).jls"))

        plot!(cc.distance,
              cc.contrast_corr,
              yscale=:log10,
              xlim=(0, NaN),
              ylim=(1e-5, 1e0),
              label="PCA $(n)",
              ylabel="5-sigma contrast",
              xlabel="radius [px]")
    end

    savefig(joinpath(paths.sequences_folder, "all_contrast_curves.png"))

end

@stage function do_adi(frames, template_psf; n_pca=[4, 8, 12, 16, 20, 24])
    paths = context["paths"]

    angle_tuple = calculate_north_angle.(frames)
    angles = map(first, angle_tuple)
    cube = framelist_to_cube(frames)


    fwhm = 2*sqrt(2*log(2))*(template_psf["SIGMA_X"] + template_psf["SIGMA_Y"]) / 2

    algs, reduced = do_pca(cube, angles, paths, n_pca)
    calc_contrast_curves(algs, cube, angles, template_psf, fwhm, paths, n_pca)
    plot_contrast_curves(paths, n_pca)

end