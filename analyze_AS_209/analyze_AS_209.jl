using Trapz

include("../src/pipeline.jl")
include("01_get_astrometry.jl")
include("02_get_photometry.jl")

function load_magnitude_data(filename)
    wl = Float64[]
    flux_density = Float64[]
    open(filename, "r") do file
        for line in eachline(file)
            if startswith(line, "#") || isempty(strip(line))
                continue
            end
            parts = split(strip(line))
            push!(wl, parse(Float64, parts[1]))
            push!(flux_density, parse(Float64, parts[2]))
        end
    end
    return wl, flux_density
end

Ks_2MASS_filename = "analyze_AS_209/2MASS_Ks.txt"
Kp_NIRC2_filename = "analyze_AS_209/NIRC2_Kp.txt"

wl_2MASS_Ks, flux_density_2MASS_Ks = load_magnitude_data(Ks_2MASS_filename)
wl_NIRC2_Kp, flux_density_NIRC2_Kp = load_magnitude_data(Kp_NIRC2_filename)
flux_density_NIRC2_Kp ./= 100

integrated_2MASS = trapz(wl_2MASS_Ks, flux_density_2MASS_Ks)  # should be close to 0.262 microns
integrated_NIRC2 = trapz(wl_NIRC2_Kp, flux_density_NIRC2_Kp)
bandpass_factor = integrated_NIRC2 / integrated_2MASS # convert from 2MASS Ks to NIRC2 Kp, assuming ideal bandpasses
@info "Integrated Bandpasses" integrated_2MASS=integrated_2MASS integrated_NIRC2=integrated_NIRC2 bandpass_factor=bandpass_factor

good_dates = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

datedict_median_coarse_locations = Dict{String, Tuple{Int,Int}}(
    "2002-06-16" => (80, 110),
    "2002-08-02" => (77, 107),
    "2002-08-21" => (77, 111),
    "2005-07-27" => (75, 117)
)

datedict_frame_coarse_locations = Dict{String, Any}(
    "2002-06-16" => ((80,111),(81,110), (79,109), (80,109)),
    "2002-08-02" => ((75,108), (76,109), (77,107),(76,108),(78,107),(78,107)),
    "2002-08-21" => ((77,111), (77, 111), (77,111), (77,111)),
    "2005-07-27" => ((76,118), (76,119), (75,118), (75,118), (76,117), (76,117), (75,118), (76,117))
)

datedict_target_key = Dict{String, String}("2002-06-16" => "AS_209_4",
                                           "2002-08-02" => "AS_209_2",
                                           "2002-08-21" => "AS_209_4",
                                           "2005-07-27" => "AS_209_2")

datedict_template_psf_keys = Dict{String, String}("2002-06-16" => "AS_209_1",
                                                  "2002-08-02" => "AS_209_1",
                                                  "2002-08-21" => "AS_209_2",
                                                  "2005-07-27" => "AS_209_1")

for date in good_dates
    @info "Processing date: $date"

    paths = ObslogPaths("/Users/jsn/landing/projects/AIR.jl/AS_209_data/", date)

    target_northup_median = load(joinpath(paths.sequences_folder, "$(datedict_target_key[date])_northup_median.fits"))
    target_northup_cube = load(joinpath(paths.sequences_folder, "$(datedict_target_key[date])_northup_cube.fits"), : )
    aligned_frames = load(joinpath(paths.sequences_folder, "$(datedict_target_key[date])_aligned_frames.fits"), :)

    if date=="2002-08-02"
        template_psf = load("AS_209_data/2002-06-16/sequences/AS_209_1_template_psf.fits")
        encircled_energy_centered_sequence = load("AS_209_data/2002-06-16/sequences/AS_209_1_centered_sequence.fits", :)
    else
        template_psf = load(joinpath(paths.sequences_folder, "$(datedict_template_psf_keys[date])_template_psf.fits"))
        encircled_energy_centered_sequence= load(joinpath(paths.sequences_folder, "$(datedict_template_psf_keys[date])_centered_sequence.fits"), :)
    end

    # high level data that everything has access to
    context = Dict{String, Any}("paths" => paths,
                                "target" => "AS_209")

    pipeline = Pipeline(save_state=false)
    add_context!(pipeline, context)

    input = [date,
            datedict_median_coarse_locations[date],
            datedict_frame_coarse_locations[date],
            target_northup_median,
            target_northup_cube,
            template_psf]
    stage_01 = Stage("01", epoch_get_astrometry; input=input)
    add_stage!(pipeline, stage_01)

    # some weirdness with one of the template PSFs due to bad flats/darks, has negative flux far away
    # so these numbers cant be too big
    r_in = 60.0
    r_out = 100.0

    kwargs = (;r_in=r_in,
               r_out=r_out)

    # correction factor using a scaling argument
    K_bandpass_2MASS = 0.262 # microns
    AS_209_band_mag = 6.961
    zero_flux = 4.283e-14 * K_bandpass_2MASS * 100*100 # W/m^2
    scaling_calibrated_flux = zero_flux * 10^(AS_209_band_mag/(-2.5)) * bandpass_factor
    @info "2MASS Calibrated Flux" calibrated_flux_W_m2=scaling_calibrated_flux

    # using color transformations from 2MASS to MKO to Kp
    J_mag_2MASS = 8.302
    H_mag_2MASS = 7.454
    Ks_mag_2MASS = 6.961

    # From Leggett et al. 2006
    K_mag_MKO = Ks_mag_2MASS + (-0.003) + (-0.025)*(J_mag_2MASS - Ks_mag_2MASS)
    H_mag_MKO = H_mag_2MASS + (-0.014) + (0.049)*(H_mag_2MASS - Ks_mag_2MASS)

    Kp_mag_MKO = K_mag_MKO + 0.22 * (H_mag_MKO - K_mag_MKO)
    @info "Kp magnitude (MKO)" Kp_mag_MKO=Kp_mag_MKO

    zeropoint_flux = 4.57E-10 * 0.351
    calibrated_flux = zeropoint_flux * 10^(Kp_mag_MKO/(-2.5))
    @info "Calibrated Flux" calibrated_flux_W_m2=calibrated_flux

    stage_02 = Stage("02", get_photometry, kwargs=kwargs)
    add_input!(stage_02, date, aligned_frames, encircled_energy_centered_sequence, calibrated_flux)
    add_stage!(pipeline, stage_02)

    run(pipeline)
end