using Printf
using OrderedCollections
using Statistics
using AstroImages
using Photometry

using AIR
import AIR.crop 

function extract_psf(northup_cropped_cube, template_psf)

    # template PSF should be normalized and have zero background
    template_psf ./= maximum(template_psf)

    initial_guess = [0.0, 0.0, 30.0, 0.0]
    search_radius = 3.0

    residuals = AstroImage[]
    extracted_psfs = AstroImage[]
    rescaled_template_psfs = AstroImage[]
    for frame in northup_cropped_cube

        derotated_template_psf = rotate_image_center(template_psf, -(frame["PARANG"]-frame["ROTPOSN"]))
        derotated_template_psf, _, _ = crop(derotated_template_psf, size(northup_cropped_cube[1]))

        residual, optimal_params = optimal_subtract_target(frame, derotated_template_psf, initial_guess, search_radius; scale_bounds=(1e-3, 1e3), offset_bounds=(-1000.0, 1000.0))

        #@info "Optimal Params" optimal_params=optimal_params

        extracted_psf = optimal_params[3] * derotated_template_psf
        rescaled_template_psf = optimal_params[3] * template_psf

        extracted_psf["FILTER"] = frame["FILTER"]

        push!(extracted_psfs, extracted_psf)
        push!(rescaled_template_psfs, rescaled_template_psf)
        push!(residuals, residual)

    end

    return extracted_psfs, rescaled_template_psfs, residuals

end

planck = 6.62607015e-34
c = 2.99792458e8

kp_lambda = 2.124e-6  # Kp filter central wavelength in meters
energy_from_wl = (wl) -> planck * c / wl
telescope_aperture_area = π * (9.96/2)^2
kp_bandpass = 0.3516  # Kp filter bandwidth in microns

function get_fluxes(psfs, template_psf)

    fluxes = Float64[]
    for psf in psfs
        aperture_radius = size(psf, 1) / 2
        aperture = CircularAperture(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, aperture_radius)
        res = photometry(aperture, psf)

        flux = energy_from_wl(kp_lambda) * res.aperture_sum/telescope_aperture_area

        push!(fluxes, flux)

    end
    @info "fluxes" fluxes=fluxes

    return fluxes

end

function get_ee_fraction(psf, r_in, r_out)

    circle_aperture = CircularAperture(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, r_out)
    annulus_aperture = CircularAnnulus(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, r_in, r_out)

    circle_photometry = photometry(circle_aperture, psf)
    annulus_photometry = photometry(annulus_aperture, psf)

    return annulus_photometry.aperture_sum / circle_photometry.aperture_sum
end

# calibrated_flux in W/m^2, has to come from something like 2MASS
@stage function get_photometry(date, aligned_frames, encircled_energy_centered_sequence, calibrated_flux, northup_cropped_cube, template_psf; save_path=nothing, r_in=60.0, r_out=100.0)

    if save_path === nothing
        save_path = joinpath("analyze_AS_209/$(date)_photometry.toml")
    end

    paths = context["paths"]

    # because some of the epochs are in PA mode, the PSF of the target rotates throughout the frame
    # must do the photometry on each individual frame rather than the median frame as a result to fit the model better

    @info "Processing date: $date"
    
    encircled_energy_fractions = Float64[]
    for psf in encircled_energy_centered_sequence
        eef = get_ee_fraction(psf, r_in, r_out)
        push!(encircled_energy_fractions, eef)
    end
    encircled_energy_fraction = mean(encircled_energy_fractions)
    encircled_energy_error = std(encircled_energy_fractions)

    annulus_fluxes = Float64[]
    for af in aligned_frames
        annulus = CircularAnnulus(size(af, 2)/2 + 0.5, size(af, 1)/2 + 0.5, r_in, r_out)
        annulus_photometry = photometry(annulus, af)
        af = annulus_photometry.aperture_sum * energy_from_wl(kp_lambda) / telescope_aperture_area / encircled_energy_fraction
        push!(annulus_fluxes, af)
    end
    @info "annulus_fluxes" annulus_fluxes=annulus_fluxes
    annulus_fluxes_mean = mean(annulus_fluxes)
    annulus_fluxes_err = std(annulus_fluxes)

    throughput_factor = calibrated_flux / annulus_fluxes_mean

    extracted_psfs, rescaled_template_psfs, residuals = extract_psf(northup_cropped_cube, template_psf)

    shot_noise = mean([sqrt(sum(psf)) for psf in rescaled_template_psfs]) * energy_from_wl(kp_lambda) / telescope_aperture_area
    background_rms(res) = sqrt(mean(abs2, parent(res)))
    background_noise = mean(background_rms.(residuals)) * energy_from_wl(kp_lambda) / telescope_aperture_area
    read_noise = get_NIRC2_readnoise(northup_cropped_cube[1]["SAMPMODE"]) * energy_from_wl(kp_lambda) / telescope_aperture_area

    fluxes = get_fluxes(rescaled_template_psfs, template_psf) # use the rescaled full versions to get all the flux in the wings that are probably suppressed in the extracted PSFs
    fluxes .*= throughput_factor

    mean_flux = mean(fluxes)
    std_flux = std(fluxes)

    ee_term = mean_flux * (encircled_energy_error/encircled_energy_fraction)
    annulus_fluxes_term = mean_flux * (annulus_fluxes_err / annulus_fluxes_mean)
    err_flux = sqrt( std_flux^2 + shot_noise^2 + background_noise^2 + read_noise^2 + ee_term^2 + annulus_fluxes_term^2 )

    spectral_flux_density = mean_flux / kp_bandpass
    err_spectral_flux_density = err_flux / kp_bandpass

    @info "Terms" encircled_energy_fraction=encircled_energy_fraction ee_term=ee_term annulus_fluxes_term=annulus_fluxes_term calibrated_flux=calibrated_flux annulus_fluxes_mean=annulus_fluxes_mean
    @info "Throughput Factor" throughput_factor=throughput_factor

    @info std_flux=std_flux shot_noise=shot_noise background_noise=background_noise read_noise=read_noise

    #@info "Final Flux" flux=mean_flux err=err_flux
    @info @sprintf("%8.3e +/- %8.3e W/m^2/um", spectral_flux_density, err_spectral_flux_density)
    @info @sprintf("%8.3e +/- %8.3e W/m^2/um", spectral_flux_density/10^(-17), err_spectral_flux_density/10^(-17)) " (1e-17 W/m^2/um units)"

    save(joinpath(paths.sequences_folder, "$(date)_extracted_psfs.fits"), extracted_psfs...)
    save(joinpath(paths.sequences_folder, "$(date)_rescaled_template_psfs.fits"), rescaled_template_psfs...)
    save(joinpath(paths.sequences_folder, "$(date)_northup_residuals.fits"), residuals...)

    save_photometry = OrderedDict{String, Any}()
    save_photometry["epoch_$(date)"] = OrderedDict(
        "flux" => mean_flux,
        "flux_units" => "W/m^2",
        "err_flux" => err_flux,
        "err_flux_units" => "W/m^2",
        "spectral_flux_density" => spectral_flux_density,
        "spectral_flux_density_units" => "W/m^2/μm",
        "err_spectral_flux_density" => err_spectral_flux_density,
        "err_spectral_flux_density_units" => "W/m^2/μm",
        "filter" => aligned_frames[1]["FILTER"],
        "filter_wl" => kp_lambda * 1e6,
        "avg_wl_units" => "μm",
        "filter_bandpass" => kp_bandpass,
        "bandpass_units" => "μm"
    )

    write_toml(save_path, save_photometry)

end