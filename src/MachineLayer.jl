"Abstract base type for machine specs"
abstract type AbstractMachineSpec end

"""
    machine_init(machine_spec, X, y)

Initialize a machine with predictor table X and target y.
"""
function machine_init end

"""
    machine_fit!(machine_spec, m)

The default implementation returns `MLJ.fit!(m, rows = m_spec.train_rows)`.
"""
function machine_fit!(m_spec::AbstractMachineSpec, m)
    return MLJ.fit!(m, rows=m_spec.train_rows)
end

"""
    machine_predict(machine_spec, m)

Produce a prediction `y_hat`.
The default implementation returns `MLJ.predict(m, X)`.
"""
function machine_predict(m_spec::AbstractMachineSpec, m, X)
    return MLJ.predict(m, X)
end

"""
    machine_complexity(machine_spec, m)

Numerical complexity of a machine.
"""
function machine_complexity end

"""
    residual_norm(machine_spec, y - y_hat)

Norm of residuals.
The default implementation returns the mean square error.
"""
function residual_norm(
    m_spec::AbstractMachineSpec,
    r::AbstractVector
    )::Real
    return dot(r, r) / length(r)
end

"""
    genome_parameter_complexity(machine_spec, p)

Numerical complexity of a parameter vector.  The default
implementation returns `machine_spec.lambda_parameter` times the
squared L2 norm of `p`.
"""
function genome_parameter_complexity(
    m_spec::AbstractMachineSpec,
    p::AbstractVector
    )::Real
    return m_spec.lambda_parameter * dot(p, p)
end

"""
    genome_complexity(machine_spec, g_spec, genome)

Numerical complexity of genome.  The default implementation
returns `machine_spec.lambda_operand` times the number of
operands in the genome.
"""
function genome_complexity(
    m_spec::AbstractMachineSpec,
    g_spec::GenomeSpec,
    genome::Genome
    )::Real
    return m_spec.lambda_operand * num_operands(genome)
end



"""
Abstract base type for solver specs, used for solving for genome parameter values.
"""
abstract type AbstractSolverSpec end

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
    parameter_solver_optimization_problem(sol_spec, g_spec, genome, optim_fn)

Build an `OptimizationProblem` around the `OptimizationFunction`.
The resulting problem object will be passed to `solve`.
The default implementation returns `OptimizationProblem(optim_fn, zeros(...), sense=MinSense)`.
"""
function parameter_solver_optimization_problem(::AbstractSolverSpec, g_spec::GenomeSpec, genome::Genome, optim_fn)
    return OptimizationProblem(optim_fn, zeros(g_spec.parameter_size), sense = MinSense)
end

"""
    parameter_solver_solve(sol_spec, optim_prob)

Solve the problem generated by `parameter_solver_optimization_problem()`.
The default implementation returns `solve(optim_prob, NelderMeade())`
"""
function parameter_solver_solve(::AbstractSolverSpec, optim_prob)
    return solve(optim_prob, NelderMeade())
end

"""
    machine_grow_and_rate(xs, y, g_spec, genome, machine_spec, sol_spec)

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
    machine_spec::AbstractMachineSpec,
    sol_spec::AbstractSolverSpec
    )::Union{Agent, Nothing}

    col_names = map(1:g_spec.output_size) do t
        "z$t"
    end

    g_c = genome_complexity(machine_spec, g_spec, genome)
    m_save = nothing

    function f(genome_parameter, _)
        last_round = run_genome(g_spec, genome, genome_parameter, xs)[end]
        outputs = map(last_round) do z
            as_vec(z, (g_spec.output_size))
        end
        Z_df = DataFrame(outputs, col_names, copycols=false)
        m = machine_init(machine_spec, Z_df, y)
        machine_fit!(machine_spec, m)
        y_hat = machine_predict(machine_spec, m)
        residuals = y - y_hat
        r_norm = residual_norm(machine_spec, residuals)
        m_c = machine_complexity(m)
        p_c = genome_parameter_complexity(machine_spec, genome_paramter)
        # This does update the outer scope `m_save`.
        m_save = (m=m, r_norm=r_norm, g_c=g_c, m_c=m_c, p_c=p_c)
        return r_norm + g_c + m_c + p_c
    end

    optim_fn = parameter_solver_optimization_function(sol_spec, f)
    optim_prob = parameter_solver_optimization_problem(sol_spec, g_spec, genome, optim_fn)
    try
        sol = parameter_solver_solve(sol_spec, optim_prob)
        if SciMLBase.successful_retcode(sol)
            @assert !isnothing(m_save)
            return Agent(r, genome, sol.u, m_save)
        else
            return nothing
        end
    catch e
        return nothing
    end
end
