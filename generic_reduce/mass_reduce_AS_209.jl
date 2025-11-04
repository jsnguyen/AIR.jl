include("../src/pipeline.jl")
include("01_sort_frames.jl")
include("02_make_masters.jl")
include("03_mass_reduce.jl")
include("04_make_table_and_plots.jl")

# know this apriori after sorting through everything
good_dates = ["2002-06-16", "2002-08-02", "2002-08-21", "2005-07-27"]

for date in good_dates
    @info "Reducing date: $date"

    paths = ObslogPaths("/Users/jsn/landing/projects/AIR.jl/AS_209_data", date)

    # high level data that everything has access to
    context = Dict{String, Any}("paths" => paths,
                                "target" => "AS_209")

    pipeline = Pipeline(save_state=false)
    add_context!(pipeline, context)

    stage_01 = Stage("01", sort_frames)
    add_stage!(pipeline, stage_01)

    stage_02 = Stage("02", make_all_masters)
    add_stage!(pipeline, stage_02)

    kwargs = (;skip_sky_sub=true,)
    stage_03 = Stage("03", mass_reduce_frames; kwargs=kwargs)
    add_stage!(pipeline, stage_03)

    stage_04 = Stage("04", make_table_and_plots)
    add_stage!(pipeline, stage_04)

    run(pipeline)

    save_context(pipeline, "generic_reduce/context/reduce_$(date)_context.jls")
end

# need to add checkpointing
# reloading specific state
# also tracking saved files and any guarantees with them...
# how to track if the outputs are written to file or not?