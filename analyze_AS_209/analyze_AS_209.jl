include("../src/pipeline.jl")
include("01_generate_sequence.jl")
include("02_make_template_psf.jl")

good_dates = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

#for date in good_dates
date = "2002-06-16"
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

run(pipeline)