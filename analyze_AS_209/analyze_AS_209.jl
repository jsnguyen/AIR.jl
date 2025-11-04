include("../src/pipeline.jl")
include("01_get_astrometry.jl")
include("02_get_photometry.jl")

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

    if date=="2002-08-02"
        template_psf = load("AS_209_data/2002-06-16/sequences/AS_209_1_template_psf_cored.fits")
    else
        template_psf = load(joinpath(paths.sequences_folder, "$(datedict_template_psf_keys[date])_template_psf.fits"))
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

    stage_02 = Stage("02", get_photometry; passthrough=true)
    add_input!(stage_02, date)
    add_stage!(pipeline, stage_02)

    run(pipeline)
end