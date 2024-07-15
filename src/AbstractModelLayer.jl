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
        agent::Agent{<:Number, <:AbstractGenome, <:AbstractVector, <:AbstractModelResult},
        xs::Vector)
    num_rows = length(xs[1])
    last_round = run_genome(g_spec, agent.genome, agent.parameter, xs)[end]
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    return model_predict(agent.extra, data_cols)
end

"""
    model_symbolic_output(g_spec, agent)

Build a symbolic form for the output of the final time step of
running `agent`'s `genome`
Then use the `agent`'s `parameter` vector,
and feed the symbolic output of the genome as input to
model result in `agent.extra` to make a prediction
in symbolic form.
"""
function model_symbolic_output(g_spec, agent)
    return model_symbolic_output(g_spec, agent.genome, agent.parameter, agent.extra)
end

"""
    model_symbolic_output(g_spec, genome, parameter, model_result)

        model_symbolic_output(g_spec, agent)

Build a symbolic form for the output of the final time step of
running `genome`
Then use the `parameter` vector,
the symbolic output of the genome,
and feed the symbolic output of the genome as input to
model result in `agent.extra` to make a prediction
in symbolic form.

Return a named tuple with lots of useful fields.

"""
function model_symbolic_output end

function model_symbolic_output(
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        parameter::AbstractVector,
        mr::AbstractModelResult)
    p, x, z = run_genome_symbolic(g_spec, genome)
    z_sym_row_mat = reshape(z, 1, :)
    y_sym = model_predict(mr, z_sym_row_mat)[1]
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
    p_subs = Dict(p[j] => parameter[j] for j in eachindex(p))
    y_sub = substitute(y_simp, p_subs)
    y_num = simplify(y_sub)

    return (p = p, x = x, z = z, p_subs = p_subs, y_sym = y_sym,
        y_lim = y_lim, y_simp = y_simp, y_sub = y_sub, y_num = y_num)
end
