include("../src/pipeline.jl")
include("01_generate_sequence.jl")
include("02_make_template_psf.jl")
include("03_adi.jl")

paths = ObslogPaths("2025-10-09", "/Users/jsn/landing/projects/AIR.jl/data/HR_8799/2025-10-09")

reduction_context = load_context("generic_reduce/reduce_2025-10-09_context.toml")
reduced_frames = load.(reduction_context["reduced_frames"])

# high level data that everything has access to
context = Dict{String, Any}("paths" => paths,
                            "target" => "HR_8799")

pipeline = Pipeline(save_state=false)
add_context!(pipeline, context)

stage_01 = Stage("01", generate_sequence_epoch, input=(reduced_frames, (449, 659) ))
add_stage!(pipeline, stage_01)

kwargs = (;coarse_size=380,
           fine_size=370,
           fixed_sigma=5.0)
sequence_index = 1

stage_02 = Stage("02", make_template_psf; kwargs=kwargs)
add_additional_input!(stage_02, (sequence_index,))
add_stage!(pipeline, stage_02)

kwargs = (;n_pca=[15, 25, 30])
stage_03 = Stage("03", do_adi, kwargs=kwargs)
add_stage!(pipeline, stage_03)

run(pipeline)