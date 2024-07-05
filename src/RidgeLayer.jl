export least_squares_ridge, least_squares_ridge_grow_and_rate
export linear_model_predict
export linear_model_symbolic_output

"""
    as_vec(u, dims)

(Private) This is a utility function to make it easier to mix
scalars and vectors in the work space when evaluating a genome.
Given a vector `u`, check that its size matches `dims` and return
it.  Given a scalar `u`, fill a vector of size `dims` with it and
return that.
"""
function as_vec end

function as_vec(u::AbstractVector, dims)
    @assert size(u)==dims "vector size is $(size(u)), expected $dims"
    return u
end

function as_vec(u::Number, dims)
    return fill(u, dims)
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
        xs::AbstractVector,
        y::AbstractVector,
        lambda::Float64,
        g_spec::GenomeSpec,
        genome::Genome,
        parameter::AbstractVector)::Tuple{Float64, AbstractVector}
    @assert g_spec.input_size == length(xs)
    outputs = run_genome(g_spec, genome, parameter, xs)
    last_round = outputs[end]
    num_output_cols = length(last_round)
    column_1s = ones(size(y, 1))
    data_cols = map(u -> as_vec(u, size(column_1s)), last_round)
    X = stack([column_1s, data_cols...])
    @assert num_output_cols + 1 == size(X, 2)
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
    least_squares_ridge_grow_and_rate(xs, y, lambda_b, lambda_p, lambda_op, g_spec, genome)

Solve for the parameter vector `p` that minimzes
`norm(y - y_hat)^2 + lambda_b * norm(b)^2 + lambda_p * norm(p)^2 + lambda_op R`,
where `y_hat` and `b` are found using `least_squares_ridge`.
`R` is the total number of operands across all instructions in `genome`.

If all goes well, return an `Agent`, whose genome is `genome`,
whose `parameter` is the best `p`, and whose `extra` is `b`.
Otherwise, return `nothing`.

"""
function least_squares_ridge_grow_and_rate(
        xs::AbstractVector,
        y::AbstractVector,
        lambda_b::Float64,
        lambda_p::Float64,
        lambda_operand::Float64,
        g_spec::GenomeSpec,
        genome::Genome)::Union{Agent, Nothing}
    u0 = zeros(g_spec.parameter_size)
    f_opt = OptimizationFunction(_LSRGR_f)
    c = _LSRGR_Context(g_spec, genome, lambda_b, xs, y, lambda_p, nothing, nothing)
    prob = OptimizationProblem(f_opt, u0, c, sense = MinSense)
    try
        sol = solve(prob, NelderMead())
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
        throw(e)
    end
end

@kwdef mutable struct _LSRGR_Context{TXs,Ty}
    g_spec::GenomeSpec
    genome::Genome
    lambda_b::Float64
    xs::TXs
    y::Ty
    lambda_p::Float64
    n::Union{Float64,Nothing}
    b::Union{Vector{Float64},Nothing}
end

function _LSRGR_f(u::Vector{Float64}, c::_LSRGR_Context{TXs,Ty}) where {TXs, Ty}
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

Return a named tuple `(p, x, w, b, y_sym, y_num)` where `p`, `x`,
and `b` are vectors of `Symbolics` objects used to represent
genome parameters, inputs, and linear model coefficients; `w` is
a vector of genome outputs in symbolic form; `y_sym` is the
symbolic representation of the linear model `dot(b, [1, w...])`;
and `y_num` is the symbolic model with subsitutions `p =
agent.parameters` and `b = agent.extra` applied, so that only
`x[j]`s remain.

"""
function linear_model_symbolic_output(
        g_spec::GenomeSpec,
        agent::Agent;
        parameter_sym = :p,
        input_sym = :x,
        coefficient_sym = :b)
    p, x, z = run_genome_symbolic(
        g_spec, agent.genome;
        parameter_sym = parameter_sym,
        input_sym = input_sym)
    b = Symbolics.variables(coefficient_sym, 1:(g_spec.output_size + 1))
    y_pred_sym = dot(b, [1, z...])
    y_pred_simp = simplify(y_pred_sym; expand = true)
    p_subs = Dict(p[j] => agent.parameter[j] for j in eachindex(p))
    b_subs = Dict(b[j] => agent.extra[j] for j in eachindex(b))
    y_num = simplify(substitute(y_pred_simp, merge(p_subs, b_subs)); expand = true)
    return (p = p, x = x, w = z, b = b, y_sym = y_pred_simp, y_num = y_num)
end

"""
    linear_model_predict(g_spec::GenomeSpec, agent::Agent, xs::Vector)

Run `agent.genome` on inputs `xs` and `agent.parameter`, and
form the linear combination of a column of 1s and the genome's
outputs using the coefficients `agent.extra`.
"""
function linear_model_predict(
        g_spec::GenomeSpec,
        agent::Agent,
        xs::Vector)
    num_rows = length(xs[1])
    last_round = run_genome(g_spec, agent.genome, agent.parameter, xs)[end]
    column_1s = ones(num_rows)
    data_cols = map(u -> as_vec(u, size(column_1s)), last_round)
    X = stack([column_1s, data_cols...])
    y_hat = X * agent.extra
    return y_hat
end
