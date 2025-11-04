using Serialization

# kwargs are for the function and are optional

# context are lightweight info like paths, target name, etc that all stages might need
# register is a way to share data between stages without passing everything through input/output, large data should be stored here or go through input/output
# context is saved in a pipeline state save, but not necessarily the register

mutable struct Stage
    name::String
    input::Vector{Any}
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
    routes::Vector{Vector{Pair{<:Union{Int, Symbol},Int}}}
end

Stage(name::String, f::Function; input::Vector=[], kwargs::NamedTuple=NamedTuple(), multithreaded::Bool=false, passthrough::Bool=false, pipeline::Union{Pipeline, Nothing}=nothing) = begin
    if pipeline !== nothing
        s = Stage(name, input, (), kwargs, (), f, Stage[], Stage[], pipeline.context, pipeline.register, multithreaded, passthrough)
        add_stage!(pipeline, s)
        return s
    else
        s = Stage(name, input, (), kwargs, (), f, Stage[], Stage[], Dict{String, Any}(), Dict{String, Any}(), multithreaded, passthrough)
        return s
    end
end

function (s::Stage)()

    if isempty(s.input) || all(isnothing, s.input)
        res = s.f(s.context, s.register; s.kwargs...)
    else
        res = s.f(s.input..., s.context, s.register; s.kwargs...)
    end

    # must be of type tuple to match struct
    if typeof(res) <: Tuple
        s.output = res
    else
        s.output = (res,)
    end

    if s.passthrough
        s.output = Tuple(s.input)
    end

    return s.output
end

Pipeline(;save_state::Bool = false) = begin
    Pipeline(Stage[], (), (), 0, Dict{Symbol, Any}(), Dict{String, Any}(), save_state, Vector{Pair{Union{Int, Symbol},Int}}())
end

function run(p::Pipeline; start_index::Int=1)

    for (i, stage) in enumerate(p.stages)

        if i < start_index
            @info "Skipping stage $(stage.name)..."
            continue
        end

        @info "Running stage $(stage.name)..."
        try
            stage()
        catch err
            io = IOContext(stderr, :limit => true, :compact => true)
            showerror(io, err, catch_backtrace())
            @error io
        end

        if !isempty(p.routes[i])
            for r in p.routes[i]
                output_index = r[1]
                dest_stage_index = r[2]

                if output_index == :all
                    add_input!(p.stages[dest_stage_index], stage.output...)
                else
                    add_input!(p.stages[dest_stage_index], stage.output[output_index])
                end

            end
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
    push!(p.routes, Pair{<:Union{Int, Symbol},Int}[])

    if length(p.stages) > 1
        connect!(p, (length(p.stages)-1, :all) => length(p.stages))
    end

    return p
end

function add_context!(p::Pipeline, context::Dict{String, Any})
    p.context = context
end

function add_input!(s::Stage, args...)
    push!(s.input, args...)
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

_context_key(expr::Symbol) = String(expr)
_context_key(expr::QuoteNode) = _context_key(expr.value)
function _context_key(expr::Expr)
    if expr.head == :.
        return _context_key(expr.args[end])          # paths.darks_file → "darks_file"
    elseif expr.head == :call
        return string(expr.args[1])                  # foo(bar) → "foo"
    end
    return sprint(show, expr)
end
_context_key(expr) = sprint(show, expr)

# Central store helper so macro stays small
function _store_in_context(ctx, key, value)
    if ctx isa AbstractDict
        ctx[key] = value
    elseif ctx === nothing
        @warn "No Dict-like `context` available to store value." key
    else
        @warn "Context is not Dict-like; skipping store." key typeof(ctx)
    end
    return value
end

# Macro: @context_store context value_expr  [or @context_store context key="foo" value_expr]
macro context_store(ctx_expr, arg1, rest...)
    key_ast = nothing
    value_expr = arg1
    trailing = rest

    # Allow optional key= override
    if arg1 isa Expr && arg1.head == :kw && arg1.args[1] == :key
        @assert !isempty(rest) "Provide a value expression after `key=`."
        key_ast = arg1.args[2]
        value_expr, trailing = rest[1], rest[2:end]
    end
    @assert isempty(trailing) "Use as: @context_store context [key=...] value_expr"

    key_expr = key_ast === nothing ? :(_context_key($(QuoteNode(value_expr)))) : esc(key_ast)

    ctx_sym = gensym(:ctx)
    value_sym = gensym(:value)

    return quote
        local $(ctx_sym) = $(esc(ctx_expr))
        local $(value_sym) = $(esc(value_expr))
        _store_in_context($(ctx_sym), $(key_expr), $(value_sym))
    end
end

"""
    connect!(pipeline, routes...)

Wire outputs from already-computed upstream stages into downstream stages

Each `route` should be a `Pair{Tuple{Int,Int},Int}` of the form
`(src_stage_idx, src_output_idx) => dest_stage_idx`.

Examples
--------
# Single connection:
connect!(pipeline, (1, 1) => 3)   # stage 1, output 1 -> stage 3

# Multiple connections in one call:
connect!(pipeline, (1, 1) => 3, (1, 2) => 2)

# Or via a collection:
routes = [(1,1)=>3, (1,2)=>2]
connect!(pipeline, routes...)
"""
function connect!(p::Pipeline, routes::Pair{<:Tuple{Int,<:Union{Int, Symbol}},Int})
    push!(p.routes[first(first(routes))], first(routes)[2] => routes[2])
end