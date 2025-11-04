include("../src/pipeline.jl")
include("01_generate_sequence.jl")
include("02_make_template_psf.jl")
include("03_adi.jl")

date = "2025-10-09"
paths = ObslogPaths("/Users/jsn/landing/projects/AIR.jl/HR_8799_data", date)

reduction_context = load_context("generic_reduce/context/reduce_2025-10-09_context.jls")
reduced_frames = load.(reduction_context["reduced_filepaths"])

# high level data that everything has access to
context = Dict{String, Any}("paths" => paths,
                            "target" => "HIP_3008")

pipeline = Pipeline(save_state=false)
add_context!(pipeline, context)

#indices = (207, 307)
#indices = (308, 417)
indices = (207, 417)

kwargs = (;rescale_frames=true)
stage_01 = Stage("01", generate_sequence_epoch, input=[reduced_frames, indices], kwargs=kwargs)
add_stage!(pipeline, stage_01)

kwargs = (;coarse_size=300,
           fine_size=200,
           fixed_sigma=5.0)
sequence_index = 1

stage_02 = Stage("02", make_template_psf; kwargs=kwargs)
add_input!(stage_02, sequence_index)
add_stage!(pipeline, stage_02)

kwargs = (;n_pca=[5, 10, 15, 20, 25, 30])
stage_03 = Stage("03", do_adi, kwargs=kwargs)
add_stage!(pipeline, stage_03)

# contrast curves arent working for some reason?


run(pipeline)