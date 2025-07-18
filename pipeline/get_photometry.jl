using Printf
using OrderedCollections
using Statistics
using AstroImages
using Photometry

using AIR
import AIR.crop 

function extract_psf(as209_northup_cropped_cube, template_psf)

    # template PSF should be normalized and have zero background
    template_psf /= maximum(template_psf)

    initial_guess = [0.0, 0.0, 30.0, 0.0]
    search_radius = 3.0

    residuals = AstroImage[]
    extracted_psfs = AstroImage[]
    rescaled_template_psfs = AstroImage[]
    for frame in as209_northup_cropped_cube

        derotated_template_psf = rotate_image_center(template_psf, -(frame["PARANG"]-frame["ROTPOSN"]))
        derotated_template_psf, _, _ = crop(derotated_template_psf, size(as209_northup_cropped_cube[1]))

        residual, optimal_params = optimal_subtract_target(frame, derotated_template_psf, initial_guess, search_radius; scale_bounds=(1e-3, 1e3), offset_bounds=(-1000.0, 1000.0))

        @info "Optimal Params" optimal_params=optimal_params

        extracted_psf = optimal_params[3] * derotated_template_psf
        rescaled_template_psf = optimal_params[3] * template_psf

        extracted_psf["FILTER"] = frame["FILTER"]
        rescaled_template_psf["FILTER"] = frame["FILTER"]

        push!(extracted_psfs, extracted_psf)
        push!(rescaled_template_psfs, rescaled_template_psf)
        push!(residuals, residual)

    end

    return extracted_psfs, rescaled_template_psfs, residuals

end

function get_fluxes(psfs)

    fluxes = Float64[]
    for psf in psfs
        aperture_radius = size(psf, 1) / 2
        aperture_area = π * aperture_radius^2
        aperture = CircularAperture(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, aperture_radius)
        res = photometry(aperture, psf)

        planck = 1.054571817e-34
        c = 2.99792458e8

        lambda = 0.0
        if occursin("Kp", psf["FILTER"])
            lambda = 2.124e-6  # Kp filter central wavelength in meters
        end

        arcsec_to_radians = 1/((1/3600)*(π/180))
        aperture_area = π*(9.96/2)^2

        average_energy = planck * c / lambda
        flux = average_energy * res.aperture_sum/aperture_area

        push!(fluxes, flux)

    end
    @info "fluxes" fluxes=fluxes

    return fluxes

end

@autolog begin

    # because some of the epochs are in PA mode, the PSF of the target rotates throughout the frame
    # must do the photometry on each individual frame rather than the median frame as a result to fit the model better

    epochs = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

    photometry = OrderedDict{String, Any}()
    for date in epochs
        @info "Processing date: $date"

        sequence_obslog = Obslog("pipeline/obslogs/$(date)_sequences.toml")

        as209_northup_cropped_cube = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_cropped_cube.fits"), :)

        # this date doesnt have a template psf, so we use the one from epoch 1
        if date == "2002-08-02"
            template_psf = load("data/2002-06-16/sequences/as209_template_psf.fits")
        else
            template_psf = load(joinpath(sequence_obslog.paths.sequences_folder, "as209_template_psf.fits"))
        end

        extracted_psfs, rescaled_template_psfs, residuals = extract_psf(as209_northup_cropped_cube, template_psf)
        fluxes = get_fluxes(rescaled_template_psfs) # use the rescaled full versions to get all the flux in the wings that are probably suppressed in the extracted PSFs

        mean_flux = mean(fluxes)
        err_flux = std(fluxes)

        @info "Final Flux" flux=mean_flux err=err_flux

        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_extracted_psfs.fits"), extracted_psfs...)
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_rescaled_template_psfs.fits"), rescaled_template_psfs...)
        save(joinpath(sequence_obslog.paths.sequences_folder, "as209_northup_residuals.fits"), residuals...)


        lambda = 0.0
        bandpass = 0.0
        if occursin("Kp", rescaled_template_psfs[1]["FILTER"])
            lambda = 2.124e-6  # Kp filter central wavelength in meters
            bandpass = 0.351e-6  # Kp filter bandwidth in meters
        end

        photometry["epoch_$(date)"] = OrderedDict(
            "flux" => mean_flux,
            "flux_units" => "W/m^2",
            "err_flux" => err_flux,
            "err_flux_units" => "W/m^2",
            "avg_wl" => lambda,
            "avg_wl_units" => "m",
            "bandpass" => bandpass,
            "bandpass_units" => "m",
            "filter" => rescaled_template_psfs[1]["FILTER"]
        )

    end

    write_toml(joinpath("pipeline/as209_photometry.toml"), photometry)

end