export AbstractMachineSpec
export machine_init, machine_fit!, machine_predict
export machine_complexity, prediction_performance, genome_parameter_complexity
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
    machine_init(mn_spec, X, y; kw_args...)

Construct a machine with predictor table X and target y.

The default implementation calls `MLJ.machine(mn_spec.machine, X, y; kw_args...)`
"""
function machine_init(mn_spec::AbstractMachineSpec, X, y; kw_args...)
    return MLJ.machine(mn_spec.model, X, y; kw_args...)
end

"""
    machine_fit!(mn_spec, m; kw_args...)

The default implementation returns `MLJ.fit!(m; kw_args...)`.
"""
function machine_fit!(mn_spec::AbstractMachineSpec, m; kw_args...)
    return MLJ.fit!(m; kw_args...)
end

"""
    machine_predict(mn_spec, m, X; kw_args...)

Produce a prediction `y_hat`.
The default implementation returns `MLJ.predict(m, X; kw_args...)`.
"""
function machine_predict(mn_spec::AbstractMachineSpec, m, X; kw_args...)
    return MLJ.predict(m, X; kw_args...)
end

"""
    machine_complexity(mn_spec, m; kw_args...)

Numerical complexity of a machine.

The default is to return 0.
"""
function machine_complexity(mn_spec::AbstractMachineSpec, m; kw_args...)
    return 0
end

"""
    prediction_performance(mn_spec, y_hat, y_ref)

Return a performance measure of the prediction `y_hat` compared to reference `y_ref`.

The default implementation applies the callable `mn_spec.performance(y_hat, y_ref)`
"""
function prediction_performance(
        mn_spec::AbstractMachineSpec,
        y_hat::AbstractVector,
        y_ref::AbstractVector
)::Real
    mn_spec.performance(y_hat, y_ref)
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
    return OptimizationProblem(
        optim_fn, zeros(g_spec.parameter_size), context, sense = MinSense)
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
    init_kw_args::Any
    fit_kw_args::Any
end

function _MGR_f(genome_parameter::AbstractVector, c::MGRContext)
    num_rows = length(c.xs[1])
    last_round = run_genome(c.g_spec, c.genome, genome_parameter, c.xs)[end]
    outputs = map(last_round) do z
        extend_if_singleton(z, num_rows)
    end
    Z_df = DataFrame(outputs, c.output_col_names, copycols = false)
    m = machine_init(c.mn_spec, Z_df, c.y; c.init_kw_args...)
    machine_fit!(c.mn_spec, m; c.fit_kw_args...)
    y_hat = machine_predict(c.mn_spec, m, Z_df)
    performance = prediction_performance(c.mn_spec, y_hat, c.y)
    m_c = machine_complexity(c.mn_spec, m)
    p_c = genome_parameter_complexity(c.mn_spec, genome_parameter)
    c.m_save = (m = m, performance = performance, m_c = m_c, p_c = p_c)
    return performance + m_c + p_c
end

@kwdef struct MachineResult{MnSpec, Mn} <: AbstractModelResult
    mn_spec::MnSpec
    m::Mn
end

"""
    machine_grow_and_rate(xs, y, g_spec, genome, mn_spec, sol_spec; init_kw_args=Dict(), fit_kw_args=())

Solve for the parameter values that produce the best machine for
predicting `y` from the outputs of the `genome` applied to `xs`.
The resulting parameter vector is stored as the `parameter` field
in the returned `Agent`.  The resulting machine is wrapped in a `MachineResult`
and stored in the `extra` field of the `Agent`.

Calls to `machine_init` will include `init_kw_args` splatted in.
Calls to `machine_fit!` will include `fit_kw_args` splatted in.

If any exception is thrown during the solving process,
the exception is suppressed, and `nothing` is returned.
"""
function machine_grow_and_rate(
        xs::AbstractVector,
        y::AbstractVector,
        g_spec::GenomeSpec,
        genome::Genome,
        mn_spec::AbstractMachineSpec,
        sol_spec::AbstractSolverSpec;
    init_kw_args::Any=NamedTuple(),
    fit_kw_args::Any=NamedTuple()
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
        nothing,
        init_kw_args,
        fit_kw_args)

    optim_fn = parameter_solver_optimization_function(sol_spec, _MGR_f)
    optim_prob = parameter_solver_optimization_problem(sol_spec, g_spec, optim_fn, c)
    try
        sol = parameter_solver_solve(sol_spec, optim_prob)
        if SciMLBase.successful_retcode(sol)
            _MGR_f(sol.u, c)
            r = sol.objective + g_c
            @assert !isnothing(c.m_save)
            return Agent(r, genome, sol.u, MachineResult(mn_spec, c.m_save.m))
        else
            return nothing
        end
    catch e
        # if isa(e, ArgumentError) || isa(e, SingularException)
        # if isa(e, SingularException)
        return nothing
        # end
        # rethrow()
    end
end

function model_predict(mr::MachineResult, input::AbstractMatrix)
    return machine_predict(mr.mn_spec, mr.m, input)
end
