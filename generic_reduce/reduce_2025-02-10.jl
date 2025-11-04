include("../src/pipeline.jl")
include("01_sort_frames.jl")
include("02_make_masters.jl")
include("03_mass_reduce.jl")
include("04_make_table_and_plots.jl")

date = "2025-02-10"
paths = ObslogPaths("/Users/jsn/landing/projects/AIR.jl/HR_8799_data", date)

# high level data that everything has access to
context = Dict{String, Any}("paths" => paths,
                            "target" => "HR_8799")

pipeline = Pipeline(save_state=false)
add_context!(pipeline, context)

stage_01 = Stage("01", sort_frames)
add_stage!(pipeline, stage_01)

stage_02 = Stage("02", make_all_masters)
add_stage!(pipeline, stage_02)

run(pipeline)

save_context(pipeline, "generic_reduce/context/reduce_$(date)_context.jls")