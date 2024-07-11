export AbstractModelResult
export model_predict
export model_symbolic_output

"""
Abstract base type for the results of fitting a model
to data.  This object can then be used for predictions.
"""
abstract type AbstractModelResult end

"""
    model_predict(mr::AbstractModelResult, X::AbstractMatrix)

Given the matrix `X` where the columns are inputs (predictors)
and the rows are points, use `mr` to predict the target value
for all the points.
"""
function model_predict end

"""
    model_predict(mr::AbstractModelResult, xs::AbstractVector{<:AbstractVector})

Stack the columns `xs` into a matrix `X` and call `model_predict`.
""" 
function model_predict(mr::AbstractModelResult, xs::AbstractVector{<:AbstractVector})
    X = stack(xs)
    return model_predict(mr, X)
end

"""
    model_predict(g_spec::GenomeSpec, agent::Agent, xs::Vector{<:AbstractVector})

Run `agent.genome` on inputs `xs` and `agent.parameter`, and
form the linear combination of the genome's
outputs using the coefficients `agent.extra`.
"""
function model_predict(
        g_spec::GenomeSpec,
        agent::Agent{<:AbstractModelResult, <:Number, <:AbstractVector, <:AbstractGenome},
        xs::Vector)
    num_rows = length(xs[1])
    last_round = run_genome(g_spec, agent.genome, agent.parameter, xs)[end]
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    return model_predict(agent.extra, data_cols)
end

"""
    model_symbolic_output(g_spec, genome; paramter_sym=:p, input_sym=:x, coefficient_sym=:b)

Build a symbolic form for the output of the final time step of
running the `genome`, and applying linear predictor coefficients.
The parameter vector, input vector, and linear predictor
coefficients are `Symbolics` objects of the form `p[j]`, `x[j]`,
and `b[j]`.  The variable names can be specified with the keyword
arguments.

Return a named tuple with lots of useful fields.

TODO Complete documentation

"""
function model_symbolic_output end
