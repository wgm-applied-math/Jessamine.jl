# This file has the code for a from-scratch implementation of
# ridge regression.  It was written so we could work on the
# genome and evolution framework before trying to interface with
# MLJ.jl.

export least_squares_ridge, least_squares_ridge_grow_and_rate

"""
    least_squares_ridge(xs, y, lambda, g_spec, genome, parameter)

Compute ridge regression.
Assume `xs` is an array of columns as predictors and `y` is a column of target values.
Apply `run_genome` using `parameter` as the parameter vector and `xs` as the inputs.
Gather a column of 1s and the output columns as a matrix `X`.
The prediction values are `y_hat = X * b`, where `b` is a column of (unknown) coefficients.
Solve for the `b` that minimizes `norm(y - y_hat)^2 + lambda * norm(b)^2`.
If all goes well, return `(norm(y - y_hat), b)`.
Otherwise return `(Inf, nothing)`.
"""
function least_squares_ridge(
        xs::AbstractArray{<:AbstractArray},
        y::AbstractArray,
        lambda::Number,
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        parameter::AbstractArray)
    @assert g_spec.input_size == length(xs)
    num_rows = length(y)
    # local outputs::Vector{Vector{Vector{Float64}}}
    outputs = run_genome(g_spec, genome, parameter, xs)
    # local last_round::Vector{Vector{Float64}}
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

Solve for the parameter vector `p` that minimizes
`norm(y - y_hat)^2 + lambda_b * norm(b)^2 + lambda_p * norm(p)^2 + lambda_op R`,
where `y_hat` and `b` are found using `least_squares_ridge`.
`R` is the total number of operands across all instructions in `genome`.

The solver starts with `p_init` for the initial value of `p`.
If `p_init` is `nothing` or not given, a vector of zeros is used.

If all goes well, return an `Agent`, whose genome is `genome`,
whose `parameter` is the best `p`, and whose `extra` is a
`BasicLinearModelResult` with coefficient vector `b`.

Otherwise, return `nothing`.

"""
function least_squares_ridge_grow_and_rate(
        xs::AbstractArray{<:AbstractArray},
        y::AbstractArray,
        lambda_b::Number,
        lambda_p::Number,
        lambda_operand::Number,
        g_spec::GenomeSpec,
        genome::AbstractGenome;
        p_init::AbstractArray = zeros(g_spec.parameter_size)
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
                return Agent(r, genome, sol.u, BasicLinearModelResult(c.b))
            end
        else
            return nothing
        end
    catch e
        if isa(e, ArgumentError) || isa(e, SingularException) || isa(e, DomainError)
            return nothing
        end
        rethrow()
    end
end

function least_squares_ridge_grow_and_rate(
        xs::AbstractArray{<:AbstractArray},
        y::AbstractArray,
        lambda_b::Number,
        lambda_p::Number,
        lambda_operand::Number,
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        p_init::Nothing
)::Union{Agent, Nothing}
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
    if isnothing(b)
        # Ridge regression failed.
        # This should be some form of infinity.
        return n
    else
        # Ridge regression succeeded.
        c.n = n
        c.b = b
        return n^2 + c.lambda_b * dot(b, b) + c.lambda_p * dot(u, u)
    end
end
