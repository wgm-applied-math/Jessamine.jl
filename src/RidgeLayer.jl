# This file has the code for a from-scratch implementation of
# ridge regression.  It was written so we could work on the
# genome and evolution framework before trying to interface with
# MLJ.jl.

export least_squares_ridge, least_squares_ridge_grow_and_rate
export linear_model_predict
export linear_model_symbolic_output

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

"""
    least_squares_ridge(xs, y, lambda, g_spec, genome, parameter)

Compute ridge regression.
Assume `xs` is an array of columns as predictors and `y` is a column of target values.
Apply `run_genome` using `parameter` as the parameter vector and `xs` as the inputs.
Gather a column of 1s and the output columns as a matrix `X`.
The prediction values are `y_hat = X * b`, where `b` is a column of (unknown) coefficients.
Solve for the `b` that minimizes `norm(y - y_hat)^2 + lambda * norm(b)^2`.
If all goes well, return `(norm(y - y_hat), b)`.
Otherwise return `nothing`.
"""
function least_squares_ridge(
        xs::Vector{Vector{Float64}},
        y::Vector{Float64},
        lambda::Float64,
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        parameter::AbstractVector)
    @assert g_spec.input_size == length(xs)
    num_rows = length(y)
    local outputs::Vector{Vector{Vector{Float64}}}
    outputs = run_genome(g_spec, genome, parameter, xs)
    local last_round::Vector{Vector{Float64}}
    last_round = outputs[end]
    num_output_cols = length(last_round)
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    X = stack(data_cols)
    @assert num_output_cols == size(X, 2)
    b = (X' * X + lambda * I) \ (X' * y)
    y_hat = X * b
    r = y - y_hat
    n = norm(r)
    if isnan(n)
        return (Inf, nothing)
    else
        return (n, b)
    end
end

"""
    least_squares_ridge_grow_and_rate(xs, y, lambda_b, lambda_p, lambda_op, g_spec, genome, p_init = zeros(...))

Solve for the parameter vector `p` that minimzes
`norm(y - y_hat)^2 + lambda_b * norm(b)^2 + lambda_p * norm(p)^2 + lambda_op R`,
where `y_hat` and `b` are found using `least_squares_ridge`.
`R` is the total number of operands across all instructions in `genome`.

The solver starts with `p_init` for the initial value of `p`.
If `p_init` is `nothing` or not given, a vector of zeros is used.   

If all goes well, return an `Agent`, whose genome is `genome`,
whose `parameter` is the best `p`, and whose `extra` is `b`.
Otherwise, return `nothing`.

"""
function least_squares_ridge_grow_and_rate(
        xs::Vector{Vector{Float64}},
        y::Vector{Float64},
        lambda_b::Float64,
        lambda_p::Float64,
        lambda_operand::Float64,
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        p_init::Vector{Float64} = zeros(g_spec.parameter_size)
)::Union{Agent, Nothing}
    optim_fn = OptimizationFunction(_LSRGR_f)
    c = _LSRGR_Context(g_spec, genome, lambda_b, xs, y, lambda_p, nothing, nothing)
    optim_prob = OptimizationProblem(optim_fn, p_init, c, sense = MinSense)
    try
        sol = solve(optim_prob, NelderMead())
        if SciMLBase.successful_retcode(sol)
            _LSRGR_f(sol.u, c)
            if isnothing(c.b)
                return nothing
            else
                r = sol.objective + lambda_operand * num_operands(genome)
                return Agent(r, genome, sol.u, c.b)
            end
        else
            return nothing
        end
    catch e
        if isa(e, ArgumentError) || isa(e, SingularException)
            return nothing
        end
        rethrow()
    end
end

function least_squares_ridge_grow_and_rate(
        xs::Vector{Vector{Float64}},
        y::Vector{Float64},
        lambda_b::Float64,
        lambda_p::Float64,
        lambda_operand::Float64,
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        p_init::Nothing
)
    return least_squares_ridge_grow_and_rate(
        xs, y, lambda_b, lambda_p, lambda_operand, g_spec, genome)
end

@kwdef mutable struct _LSRGR_Context{TXs, Ty}
    g_spec::GenomeSpec
    genome::AbstractGenome
    lambda_b::Float64
    xs::TXs
    y::Ty
    lambda_p::Float64
    n::Union{Float64, Nothing}
    b::Union{Vector{Float64}, Nothing}
end

function _LSRGR_f(u::Vector{Float64}, c::_LSRGR_Context{TXs, Ty}) where {TXs, Ty}
    n, b = least_squares_ridge(c.xs, c.y, c.lambda_b, c.g_spec, c.genome, u)
    c.n = n
    c.b = b
    return n^2 + c.lambda_b * dot(b, b) + c.lambda_p * dot(u, u)
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
        agent::Agent{<:Number, <:AbstractVector, <:AbstractVector, <:AbstractGenome};
        parameter_sym = :p,
        input_sym = :x,
        coefficient_sym = :b)
    p, x, z = run_genome_symbolic(
        g_spec, agent.genome;
        parameter_sym = parameter_sym,
        input_sym = input_sym)
    b = Symbolics.variables(coefficient_sym, 1:(g_spec.output_size))
    y_sym = dot(z, b)
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
    b_subs = Dict(b[j] => agent.extra[j] for j in eachindex(b))
    y_sub = substitute(y_simp, merge(p_subs, b_subs))
    y_num = simplify(y_sub)
    return (p = p, x = x, z = z, b = b, p_subs = p_subs, b_subs = b_subs, y_sym = y_sym,
        y_lim = y_lim, y_simp = y_simp, y_sub = y_sub, y_num = y_num)
end

"""
    linear_model_predict(g_spec::GenomeSpec, agent::Agent, xs::Vector)

Run `agent.genome` on inputs `xs` and `agent.parameter`, and
form the linear combination of the genome's
outputs using the coefficients `agent.extra`.
"""
function linear_model_predict(
        g_spec::GenomeSpec,
        agent::Agent{<:Number, <:AbstractVector, <:AbstractVector, <:AbstractGenome},
        xs::Vector)
    num_rows = length(xs[1])
    last_round = run_genome(g_spec, agent.genome, agent.parameter, xs)[end]
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    X = stack(data_cols)
    y_hat = X * agent.extra
    return y_hat
end
