
# kwargs are for the function and are optional
# additional_input is some kind of manual input that the user provides that are necessary

# context are lightweight info like paths, target name, etc that all stages might need
# register is a way to share data between stages without passing everything through input/output, large data should be stored here or go through input/output
# context is saved in a pipeline state save, but not necessarily the register

mutable struct Stage
    name::String
    input::Tuple
    additional_input::Tuple
    output::Tuple
    kwargs::NamedTuple
    byproducts::Tuple
    f::Function
    next::Vector{Stage}
    prev::Vector{Stage}
    context::Any
    register::Dict{String, Any}
    multithreaded::Bool
    passthrough::Bool
end

mutable struct Pipeline
    stages::Vector{Stage}
    final_output::Tuple
    last_output::Tuple
    last_completed::Int
    context::Any
    register::Dict{String, Any}
    save_state::Bool
end

Stage(name::String, f::Function; input::Tuple=(), additional_input::Tuple=(), kwargs::NamedTuple=NamedTuple(), multithreaded::Bool=false, passthrough::Bool=false, pipeline::Union{Pipeline, Nothing}=nothing) = begin
    if pipeline !== nothing
        s = Stage(name, input, additional_input, (), kwargs, (), f, Stage[], Stage[], pipeline.context, pipeline.register, multithreaded, passthrough)
        add_stage!(pipeline, s)
        return s
    else
        return Stage(name, input, additional_input, (), kwargs, (), f, Stage[], Stage[], Dict{String, Any}(), Dict{String, Any}(), multithreaded, passthrough)
    end
end

function (s::Stage)()

    # can add additional input
    full_input = (s.input..., s.additional_input...)

    if length(full_input) == 0 || all(isnothing, full_input)
        res = s.f(s.context, s.register; s.kwargs...)
    else
        res = s.f(full_input..., s.context, s.register; s.kwargs...)
    end

    if typeof(res) <: Tuple
        s.output = res
    else
        s.output = (res,)
    end

    if s.passthrough
        s.output = s.input
    end

    return s.output
end

Pipeline(;save_state::Bool = false) = begin
    Pipeline(Stage[], (), (), 0, Dict{Symbol, Any}(), Dict{String, Any}(), save_state)
end

function run(p::Pipeline; start_index::Int=1)

    for (i, stage) in enumerate(p.stages)

        if i < start_index
            @info "Skipping stage $(stage.name)..."
            continue
        end

        if i != 1
            stage.input = p.stages[i-1].output
        end
        
        @info "Running stage $(stage.name)..."
        try
            stage()
        catch e
            @error "Error in stage $(stage.name): $e"
            throw(e)
        end

        p.last_completed = i
        p.last_output = stage.output
        p.register[stage.name] = stage.output

        if p.save_state
            save_state(p, "pipeline_state.jls")
        end

    end
end

function add_stage!(p::Pipeline, s::Stage)
    s.context = p.context
    s.register = p.register
    push!(p.stages, s)
    return p
end

function add_context!(p::Pipeline, context::Dict{String, Any})
    p.context = context
end

function add_additional_input!(s::Stage, input::Tuple)
    s.additional_input = input
end

# serializing the whole thing can mean huge files...
function save_state(p::Pipeline, filepath::String)
    open(filepath, "w") do io
        serialize(io, p)
    end
end

function load_state(filepath::String)
    open(filepath, "r") do io
        return deserialize(io)
    end
end

function save_context(p::Pipeline, filepath::String)
    open(filepath, "w") do io
        serialize(io, p.context)
    end
end

function load_context(filepath::String)
    open(filepath, "r") do io
        return deserialize(io)
    end
end

macro stage(def)
    @assert def.head == :function

    sig, body = def.args

    ctx = gensym(:context)
    reg = gensym(:register)
    fn = sig.args[1]
    args = sig.args[2:end]

    new_sig = Expr(:call, fn, args..., ctx, reg)

    new_body = Expr(:block,
                    :(local context = $ctx),
                    :(local register = $reg),
                    body)

    return esc(Expr(:function, new_sig, new_body))
end

"""
    @context_save f(dest, args...)

Wrap a save-like call and append `string(dest)` into `context["saved_files"]`.
Requires being called inside a `@stage`-wrapped function so `context` is in scope.
"""
macro context_save(arg1, rest...)
    key_expr = :( "saved_files" )
    call_expr = arg1
    trailing = rest

    if arg1 isa Expr && arg1.args[1] == :key && (arg1.head == :kw || arg1.head === :(=))
        @assert !isempty(rest) "Provide a save call after `key=`."
        key_expr = arg1.args[2]
        call_expr, trailing = rest[1], rest[2:end]
    end

    @assert isempty(trailing) "Use as: @context_save [key=...] f(dest, ...)"
    @assert call_expr isa Expr && call_expr.head === :call "Use as: @context_save f(dest, ...)"

    fn = call_expr.args[1]
    args = call_expr.args[2:end]
    @assert !isempty(args) "First argument must be the destination (e.g., path/URI)."

    dest_expr  = args[1]
    other_args = args[2:end]

    dest_sym = gensym(:dest)
    res_sym  = gensym(:res)
    key_sym  = gensym(:key)

    call_ast = Expr(:call, esc(fn), dest_sym, map(esc, other_args)...)

    return quote
        local $(key_sym)  = $(esc(key_expr))
        local $(dest_sym) = $(esc(dest_expr))
        local $(res_sym)  = $(call_ast)

        try
            if isdefined(Main, :context) && context !== nothing && context isa AbstractDict
                local _list = get!(context, $(key_sym), String[])
                push!(_list, string($(dest_sym)))
                context[$(key_sym)] = _list
            else
                @warn "No Dict-like `context` available to log saved file."
            end
        catch _e
            @warn "Could not update context[$(key_sym)]: $_e"
        end

        $(res_sym)
    end
end


# add serialization of input/output later
# should be able to save pipeline state
