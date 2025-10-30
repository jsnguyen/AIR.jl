include("../src/pipeline.jl")
include("01_sort_frames.jl")
include("02_make_masters.jl")
include("03_mass_reduce.jl")
include("04_make_table_and_plots.jl")

paths = ObslogPaths("2025-10-09", "/Users/jsn/landing/projects/AIR.jl/HR_8799_data/2025-10-09")

# high level data that everything has access to
context = Dict{String, Any}("paths" => paths,
                            "target" => "HR_8799")

pipeline = Pipeline(save_state=false)
add_context!(pipeline, context)

stage_01 = Stage("01", sort_frames)
add_stage!(pipeline, stage_01)

stage_02 = Stage("02", make_all_masters)
add_stage!(pipeline, stage_02)

# fix this
kwargs = (;guest_obslog = Obslog("generic_reduce/obslogs/01_2025-02-10_obslog.toml"),
          use_guest_flats = true,
          use_guest_skies = true)

stage_03 = Stage("03", mass_reduce_frames; kwargs=kwargs)
add_stage!(pipeline, stage_03)

stage_04 = Stage("04", make_table_and_plots)
add_stage!(pipeline, stage_04)

run(pipeline)

# need to add checkpointing
# reloading specific state
# also tracking saved files and any guarantees with them...
# how to track if the outputs are written to file or not?