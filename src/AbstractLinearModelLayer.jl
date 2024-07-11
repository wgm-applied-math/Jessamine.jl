export AbstractLinearModelResult
export BasicLinearModelResult, coefficients, intercept
export linear_model_predict
export linear_model_symbolic_output

export extend_if_singleton

"""
Abstract base type for the results of fitting linear models.
"""
abstract type AbstractLinearModelResult end

"""
A `BasicLinearModelResult` has a vector of coefficients
plust a bias or intercept.

"""
@kwdef struct BasicLinearModelResult <: AbstractLinearModelResult
    coefficients::Vector{Float64}
    intercept::Float64
end

BasicLinearModelResult(b) = BasicLinearModelResult(b, 0.0)

function coefficients(r::BasicLinearModelResult)
    return r.coefficients
end

function intercept(r::BasicLinearModelResult)
    return r.intercept
end

"""
    linear_model_symbolic_output(g_spec, genome; paramter_sym=:p, input_sym=:x, coefficient_sym=:b)

Build a symbolic form for the output of the final time step of
running the `genome`, and applying linear predictor coefficients.
The parameter vector, input vector, and linear predictor
coefficients are `Symbolics` objects of the form `p[j]`, `x[j]`,
and `b[j]`.  The variable names can be specified with the keyword
arguments.

Return a named tuple with lots of useful fields.

TODO Complete documentation

"""
function linear_model_symbolic_output(
        g_spec::GenomeSpec,
        agent::Agent{<:AbstractLinearModelResult, <:Number, <:AbstractVector, <:AbstractGenome};
        parameter_sym = :p,
        input_sym = :x,
        coefficient_sym = :b)
    p, x, z = run_genome_symbolic(
        g_spec, agent.genome;
        parameter_sym = parameter_sym,
        input_sym = input_sym)
    # This includes `b[1] =`` $b_0$ as an intercept or bias term.
    # Which throws off the numbering.
    b = Symbolics.variables(coefficient_sym, 0:(g_spec.output_size))
    y_sym = dot(z, b[2:end]) + b[1]
    used_vars = Set(v.name for v in Symbolics.get_variables(y_sym))
    # To handle rational functions that have things like 1/(x/0),
    # replace Inf with W and do a limit as W -> Inf.
    # First, grind through and make sure we have a unique symbol.
    j = 0
    local W
    while true
        W = Symbolics.variable(:W, j)
        if !(W in used_vars)
            break
        end
        j += 1
    end
    y_W = substitute(y_sym, Dict([Inf => W]))
    y_lim = Symbolics.limit(y_W.val, W.val, Inf)
    y_simp = simplify(y_lim)
    p_subs = Dict(p[j] => agent.parameter[j] for j in eachindex(p))
    b_num = coefficients(agent.extra)
    b_subs = Dict(b[1+j] => b_num[j] for j in eachindex(b_num))
    # This is really $b_0$:
    b_subs[b[1]] = intercept(agent.extra)
    y_sub = substitute(y_simp, merge(p_subs, b_subs))
    y_num = simplify(y_sub)
    return (p = p, x = x, z = z, b = b, p_subs = p_subs, b_subs = b_subs, y_sym = y_sym,
        y_lim = y_lim, y_simp = y_simp, y_sub = y_sub, y_num = y_num)
end

"""
    linear_model_predict(lmr::AbstractLinearModelResult, X::Matrix)

Return `X * coefficients(lmr) + intercept(lmr)`
"""
function linear_model_predict(lmr::AbstractLinearModelResult, X::Matrix)
    return X * coefficients(lmr) .+ intercept(lmr)
end

"""
    linear_model_predict(lmr::AbstractLinearModelResult, xs::AbstractVector{<:AbstractVector})

Return stack the columns `xs` into a matrix `X` and call `linear_model_predict`.
""" 
function linear_model_predict(lmr::AbstractLinearModelResult, xs::AbstractVector{<:AbstractVector})
    X = stack(xs)
    return linear_model_predict(lmr, X)
end

"""
    linear_model_predict(g_spec::GenomeSpec, agent::Agent, xs::Vector{<:AbstractVector})

Run `agent.genome` on inputs `xs` and `agent.parameter`, and
form the linear combination of the genome's
outputs using the coefficients `agent.extra`.
"""
function linear_model_predict(
        g_spec::GenomeSpec,
        agent::Agent{<:AbstractLinearModelResult, <:Number, <:AbstractVector, <:AbstractGenome},
        xs::Vector)
    num_rows = length(xs[1])
    last_round = run_genome(g_spec, agent.genome, agent.parameter, xs)[end]
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    return linear_model_predict(agent.extra, data_cols)
end

"""
    extend_if_singleton(v::AbstractVector, m::Int)

If `v` is a singleton, as in `v = [x]`,
return `[x, x, ...]`, that is, a vector of `m` copies of `x`.
Otherwise, assert that `length(v) == m` and return `v`.
So the result is always a vector of length `m`.
"""
function extend_if_singleton(v::AbstractVector, m::Int)
    if length(v) == 1
        return fill(v[1], m)
    else
        @assert length(v)==m "m = $m, v = $v"
        return v
    end
end

function extend_if_singleton(v::AbstractVector, shape::Tuple{Int})
    (m,) = shape
    return extend_if_singleton(v, m)
end
