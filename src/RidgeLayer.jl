export least_squares_ridge, least_squares_ridge_grow_and_rate
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
    @assert size(u) == dims
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
    offset_col = ones(size(y, 1))
    data_cols = map(u -> as_vec(u, size(offset_col)), last_round)
    X = stack([offset_col, data_cols...])
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

    function f(u, _)
        n, b = least_squares_ridge(xs, y, lambda_b, g_spec, genome, u)
        return n^2 + lambda_b * norm(b)^2 + lambda_p * norm(u)^2
    end

    f_opt = OptimizationFunction(f)
    prob = OptimizationProblem(f_opt, u0, sense = MinSense)
    sol = solve(prob, NelderMead())
    if SciMLBase.successful_retcode(sol)
        n, b = least_squares_ridge(xs, y, lambda_b, g_spec, genome, sol.u)
        if isnothing(b)
            return nothing
        else
            r = sol.objective + lambda_operand * num_operands(genome)
            return Agent(r, genome, sol.u, b)
        end
    else
        return nothing
    end
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
        agent::Agent{<:AbstractVector};
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
