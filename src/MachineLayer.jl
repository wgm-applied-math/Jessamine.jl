"Abstract base type for machine specs"
abstract type AbstractMachineSpec end

"""
    machine_init(machine_spec, X, y)

Initialize a machine with predictor table X and target y."""
function machine_init end

"""
    machine_fit!(machine_spec, m)

Fit or otherwise learn."""
function machine_fit! end

"""
    machine_predict(machine_spec, m)

Produce a prediction `y_hat`.
"""
function machine_predict end

"""
    residual_norm(machine_spec, y - y_hat)

Norm of residuals."""
function residual_norm end

"""
    machine_complexity(machine_spec, m)

Numerical complexity of a machine.
"""
function machine_complexity end

"""
    genome_complexity(machine_spec, g_spec, genome)

Numerical complexity of genome.
"""
function genome_complexity end

"""
    genome_parameter_complexity(machine_spec, p)

Numerical complexity of a parameter vector.
"""
function genome_parameter_complexity end

function machine_grow_and_rate(
    xs::AbstractVector,
    y::AbstractVector,
    lambda_c::Real,
    g_spec::GenomeSpec,
    genome::Genome,
    machine_spec::AbstractMachineSpec
    ;
    genome_parameter_init=zeros(g_spec.parameter_size),
    genome_parameter_solver=NelderMead(),
    )::Union{Agent, Nothing}

    col_names = map(1:g_spec.output_size) do t
        "z$t"
    end

    g_c = genome_complexity(machine_spec, g_spec, genome)
    m_save = [nothing]

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
        m_save[1] = (m=m, r_norm=r_norm, g_c=g_c, m_c=m_c, p_c=p_c)
        return r_norm + g_c + m_c + p_c
    end

    f_opt = OptimizationFunction(f)
    prob = OptimizationProblem(f_opt, genome_parameter_init, sense = MinSense)
    try
        sol = solve(prob, parameter_solver)
        if SciMLBase.successful_retcode(sol)
            return Agent(r, genome, sol.u, m_save[1])
        else
            return nothing
        end
    catch e
        return nothing
    end
end
