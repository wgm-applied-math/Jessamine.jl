export AbstractMachineSpec
export machine_init, machine_fit!, machine_predict
export machine_complexity, residual_norm, genome_parameter_complexity
export genome_complexity
export AbstractSolverSpec, DefaultSolverSpec
export parameter_solver_optimization_function
export parameter_solver_optimization_problem
export parameter_solver_solve
export machine_grow_and_rate
export MGRContext


"Abstract base type for machine specs"
abstract type AbstractMachineSpec end

"""
    machine_init(mn_spec, X, y)

Construct a machine with predictor table X and target y.
"""
function machine_init end

"""
    machine_fit!(mn_spec, m)

The default implementation returns `MLJ.fit!(m, rows = m_spec.train_rows)`.
"""
function machine_fit!(mn_spec::AbstractMachineSpec, m)
    return MLJ.fit!(m, rows = mn_spec.train_rows)
end

"""
    machine_predict(mn_spec, m, X)

Produce a prediction `y_hat`.
The default implementation returns `MLJ.predict(m, X)`.
"""
function machine_predict(mn_spec::AbstractMachineSpec, m, X)
    return MLJ.predict(m, X)
end

"""
    machine_complexity(mn_spec, m)

Numerical complexity of a machine.
"""
function machine_complexity end

"""
    residual_norm(mn_spec, y - y_hat)

Norm of residuals.
The default implementation returns the mean square error.
"""
function residual_norm(
        mn_spec::AbstractMachineSpec,
        r::AbstractVector
)::Real
    return dot(r, r) / length(r)
end

"""
    genome_parameter_complexity(mn_spec, p)

Numerical complexity of a parameter vector.  The default
implementation returns `mn_spec.lambda_parameter` times the
squared L2 norm of `p`.
"""
function genome_parameter_complexity(
        mn_spec::AbstractMachineSpec,
        p::AbstractVector
)::Real
    return mn_spec.lambda_parameter * dot(p, p)
end

"""
    genome_complexity(mn_spec, g_spec, genome)

Numerical complexity of genome.  The default implementation
returns `mn_spec.lambda_operand` times the number of
operands in the genome.
"""
function genome_complexity(
        mn_spec::AbstractMachineSpec,
        g_spec::GenomeSpec,
        genome::Genome
)::Real
    return mn_spec.lambda_operand * num_operands(genome)
end

"""
Abstract base type for solver specs, used for solving for genome parameter values.
"""
abstract type AbstractSolverSpec end

"""
Use this if you want all default solving procedures.
"""
struct DefaultSolverSpec <: AbstractSolverSpec
end


"""
    parameter_solver_optimization_function(sol_spec, f)

Build an `OptimizationFunction` around the function `f(genome_parameter_vector, _)`.
The resulting objective function object will eventually be passed to `solve()`.
The default implementation returns `OptimizationFunction(f)`.
"""
function parameter_solver_optimization_function(::AbstractSolverSpec, f)
    return OptimizationFunction(f)
end

"""
    parameter_solver_optimization_problem(sol_spec, g_spec, optim_fn, context)

Build an `OptimizationProblem` around the `OptimizationFunction`.
The resulting problem object will be passed to `solve`.
The default implementation returns `OptimizationProblem(optim_fn, zeros(...), context, sense=MinSense)`.
"""
function parameter_solver_optimization_problem(
        ::AbstractSolverSpec, g_spec::GenomeSpec, optim_fn, context)
    return OptimizationProblem(optim_fn, zeros(g_spec.parameter_size), context, sense = MinSense)
end

"""
    parameter_solver_solve(sol_spec, optim_prob)

Solve the problem generated by `parameter_solver_optimization_problem()`.
The default implementation returns `solve(optim_prob, NelderMeade())`
"""
function parameter_solver_solve(::AbstractSolverSpec, optim_prob)
    return solve(optim_prob, NelderMead())
end


@kwdef mutable struct MGRContext{MnSpec, TXs, Ty}
    g_spec::GenomeSpec
    genome::AbstractGenome
    mn_spec::MnSpec
    xs::TXs
    y::Ty
    output_col_names::Vector{String}
    m_save::Any
end

function _MGR_f(genome_parameter::AbstractVector, c::MGRContext)
    num_rows = length(c.xs[1])
    last_round = run_genome(c.g_spec, c.genome, genome_parameter, c.xs)[end]
    outputs = map(last_round) do z
        extend_if_singleton(z, num_rows)
    end
    Z_df = DataFrame(outputs, c.output_col_names, copycols = false)
    m = machine_init(c.mn_spec, Z_df, c.y)
    machine_fit!(c.mn_spec, m)
    y_hat = machine_predict(c.mn_spec, m, Z_df)
    residuals = c.y - y_hat
    r_norm = residual_norm(c.mn_spec, residuals)
    m_c = machine_complexity(c.mn_spec, m)
    p_c = genome_parameter_complexity(c.mn_spec, genome_parameter)
    c.m_save = (m = m, r_norm = r_norm, m_c = m_c, p_c = p_c)
    return r_norm + m_c + p_c
end

"""
    machine_grow_and_rate(xs, y, g_spec, genome, mn_spec, sol_spec)

Solve for the parameter values that produce the best machine for
predicting `y` from the outputs of the `genome` applied to `xs`.
The resulting parameter vector is stored as the `parameter` field
in the returned `Agent`.  A named typle is stored as the `extra`
field in the returned `Agent` with these fields:

- `m`: the machine
- `r_norm`: the norm of the residuals from using `m` to predict `y`
- `g_c`: complexity of the `genome`
- `m_c`: complexity of the machine `m`
- `p_c`: complexity of the parameter vector
"""
function machine_grow_and_rate(
        xs::AbstractVector,
        y::AbstractVector,
        g_spec::GenomeSpec,
        genome::Genome,
        mn_spec::AbstractMachineSpec,
        sol_spec::AbstractSolverSpec
)::Union{Agent, Nothing}
    output_col_names = map(1:(g_spec.output_size)) do t
        "z$t"
    end

    g_c = genome_complexity(mn_spec, g_spec, genome)

    c = MGRContext(
        g_spec,
        genome,
        mn_spec,
        xs,
        y,
        output_col_names,
        nothing)

    optim_fn = parameter_solver_optimization_function(sol_spec, _MGR_f)
    optim_prob = parameter_solver_optimization_problem(sol_spec, g_spec, optim_fn, c)
    try
        sol = parameter_solver_solve(sol_spec, optim_prob)
        if SciMLBase.successful_retcode(sol)
            _MGR_f(sol.u, c)
            r = sol.objective + g_c
            @assert !isnothing(c.m_save)
            return Agent(r, genome, sol.u, (c.m_save ..., g_c = g_c))
        else
            return nothing
        end
    catch e
        # if isa(e, ArgumentError) || isa(e, SingularException)
        if isa(e, SingularException)
            return nothing
        end
        rethrow()
    end
end
