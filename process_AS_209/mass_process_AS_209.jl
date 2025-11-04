include("../src/pipeline.jl")
include("01_generate_sequence.jl")
include("02_make_template_psf.jl")
include("03_register_frames.jl")
include("04_rdi_subtract.jl")

good_dates = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

datedict_target_indices = Dict{String, Dict{String, Tuple{Int, Int}}}(
                          "2002-06-16" => Dict{String, Tuple{Int, Int}}("HBC_630"  => (412, 440),
                                                                        "HBC_650"  => (441, 465),
                                                                        "AS_209"   => (466, 490)),
                          "2002-08-02" => Dict{String, Tuple{Int, Int}}("AS_209"   => (78 , 109)),
                          "2002-08-21" => Dict{String, Tuple{Int, Int}}("AS_209"   => (84 , 117),
                                                                        "TYC_2307" => (192, 219)),
                          "2005-07-27" => Dict{String, Tuple{Int, Int}}("AS_209"   => (181, 205),
                                                                        "T222007"  => (373, 398))
)

datedict_unsaturated_sequences = Dict{String, Tuple}("2002-06-16" => ("AS_209_1", "HBC_630_1"),
                                                     "2002-08-02" => (),
                                                     "2002-08-21" => ("AS_209_2", "TYC_2307_1"),
                                                     "2005-07-27" => ("AS_209_1", "T222007_1"))

datedict_template_psf_keys = Dict{String, String}("2002-06-16" => "AS_209_1",
                                                  "2002-08-02" => "AS_209_1",
                                                  "2002-08-21" => "AS_209_2",
                                                  "2005-07-27" => "AS_209_1")

datedict_target_keys = Dict{String, Tuple}("2002-06-16" => ("AS_209_4", "HBC_630_3"),
                                           "2002-08-02" => ("AS_209_2", "AS_209_2"),
                                           "2002-08-21" => ("AS_209_4", "TYC_2307_3"),
                                           "2005-07-27" => ("AS_209_2", "T222007_4"))

datedict_sizes = Dict{String, Tuple}("2002-06-16" => (530, 500, 400),
                                     "2002-08-02" => (450, 430, 400),
                                     "2002-08-21" => (420, 410, 400),
                                     "2005-07-27" => (420, 410, 400))

for date in good_dates
    @info "Processing date: $date"

    paths = ObslogPaths("/Users/jsn/landing/projects/AIR.jl/AS_209_data/", date)

    reduction_context = load_context("generic_reduce/context/reduce_$(date)_context.jls")
    reduced_frames = load.(reduction_context["reduced_filepaths"])

    # high level data that everything has access to
    context = Dict{String, Any}("paths" => paths,
                                "target" => "AS_209")

    pipeline = Pipeline(save_state=false)
    add_context!(pipeline, context)

    stage_01 = Stage("01", generate_sequence_epoch, input=[datedict_target_indices[date], reduced_frames])
    add_stage!(pipeline, stage_01)
    connect!(pipeline, (1, 1) => 3)

    kwargs = (;coarse_size=300,
            fine_size=256,
            fixed_sigma=2.0,
            quantile_threshold=0.99999,
            use_cored=true)
    stage_02 = Stage("02", make_template_psf; kwargs=kwargs)
    add_input!(stage_02, datedict_unsaturated_sequences[date])
    add_stage!(pipeline, stage_02)


    injected_template_psf = nothing
    if date=="2002-08-02"
        injected_template_psf = load("AS_209_data/2002-06-16/sequences/AS_209_1_template_psf_cored.fits")
    end

    kwargs = (;sizes=datedict_sizes[date],
            injected_template_psf=injected_template_psf)
    stage_03 = Stage("03", register_sequences; kwargs=kwargs)
    add_input!(stage_03, datedict_template_psf_keys[date], datedict_target_keys[date])
    add_stage!(pipeline, stage_03)

    stage_04 = Stage("04", rdi_subtract)
    add_input!(stage_04, datedict_target_keys[date])
    add_stage!(pipeline, stage_04)

    try
        run(pipeline)
    catch err
        @info pipeline.routes
        io = IOContext(stderr, :limit => true, :compact => true)
        showerror(io, err, catch_backtrace())
        rethrow(err)
    end
end