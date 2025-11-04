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

function do_pca(cube, angles, paths, target, n_pca)

    algs = ADI.ADIAlgorithm[]
    reduced = AstroImage[]

    Threads.@threads for i in eachindex(n_pca)
        n = n_pca[i]
        @info "Calculating PCA reduction for PCA n = $(n)..."

        alg = PCA(n)
        res = AstroImage(alg(cube, angles))

        res["N_PCA"] = n

        push!(algs, alg)
        push!(reduced, res)

        save(joinpath(paths.sequences_folder, "$(target)_reduced_pca_$(n).fits"), res)

    end

    for i in eachindex(n_pca)
        xs, ys = axes(reduced[i])
        heatmap(xs, ys, transpose(reduced[i]); aspect_ratio=1, xlim=extrema(xs), ylim=extrema(ys), clims=(-20,20))
        savefig(joinpath(paths.sequences_folder, "$(target)_reduced_pca_$(n_pca[i]).png"))
    end

    return algs, reduced

end

function calc_contrast_curves(algs, cube, angles, template_psf, fwhm, paths, target, n_pca)

    Threads.@threads for i in eachindex(n_pca)
        n = n_pca[i]
        @info "Calculating contrast curve for PCA n = $(n)..."
        cc = contrast_curve(algs[i], cube, angles, template_psf; fwhm=fwhm, snr=5.0)
        serialize(joinpath(paths.sequences_folder, "$(target)_contrast_curve_pca_$(n).jls"), cc)
    end

end

function plot_contrast_curves(target, paths, n_pca)

    plot()

    for n in n_pca
        cc = deserialize(joinpath(paths.sequences_folder, "$(target)_contrast_curve_pca_$(n).jls"))

        plot!(cc.distance,
              cc.contrast_corr,
              yscale=:log10,
              xlim=(0, NaN),
              ylim=(1e-8, 1e0),
              label="PCA $(n)",
              ylabel="5-sigma contrast",
              xlabel="radius [px]")
        title!("Contrast Curves NPCA=$(n)")
    end

    savefig(joinpath(paths.sequences_folder, "$(target)_all_contrast_curves.png"))

end

@stage function do_adi(frames, template_psf; n_pca=[4, 8, 12, 16, 20, 24])
    paths = context["paths"]

    angle_tuple = calculate_north_angle.(frames)
    angles = map(first, angle_tuple)
    cube = framelist_to_cube(frames)

    fwhm = 2*sqrt(2*log(2))*(template_psf["SIGMA_X"] + template_psf["SIGMA_Y"]) / 2

    algs, reduced = do_pca(cube, angles, paths, context["target"], n_pca)
    calc_contrast_curves(algs, cube, angles, template_psf, fwhm, paths, context["target"], n_pca)
    plot_contrast_curves(context["target"], paths, n_pca)

end