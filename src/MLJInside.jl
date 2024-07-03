# For using MLJ machinery applied to the output of a genome

using MLJLinearModels

# Some general defaults...

"""
    residual_norm(m_spec, r)

In general, return the mean square error.
"""
function residual_norm(m_spec::AbstractMachineSpec, r::AbstractVector)
    return sum(r .* r) / length(r)
end

"""
    genome_complexity(m_spec, g_spec, genome)

In general, return `m_spec.lambda_operand` times the number of
operands in the genome.
"""

function genome_complexity(
    m_spec::AbstractMachineSpec,
    g_spec::GenomeSpec,
    genome::Genome)
    return m_spec.lambda_operand * num_operands(genome)
end

"""
    genome_parameter_complexity(m_spec, p)

In general, return `m_spec.lambda_parameter` times the L2 norm of `p`.
"""
function genome_parameter_complexity(m_spec::AbstractMachineSpec, p::AbstractVector)
    return
end


"""
    machine_fit!(m_spec, m)

In general, return `MLJ.fit!(m, rows = m_spec.train_rows)`.
"""
function machine_fit!(m_spec::AbstractMachineSpec, m)
    return MJL.fit!(m, rows=m_spec.train_rows)
end

"""
    machine_predict(m_spec, mach, X)

In general, return `MLJ.predict(m, X)`.
"""

function machine_predict(m_spec::AbstractMachineSpec, m, X)
    return MLJ.predict(m, X)
end


@kwdef struct RidgeRegressionMachineSpec <: AbstractMachineSpec
    include_bias::Bool
    lambda_ridge::Real
    lambda_parameter::Real
    lambda_operand::Real
    train_rows::AbstractVector
    test_rows::AbstractVector
end

function machine_init(m_spec::RidgeRegressionMachineSpec, X, y)
    return MJL.machine(
        RidgeRegressor(
            lambda=m_spec.lambda_ridge,
            bias=m_spec.include_bias),
        X, y)
end

function machine_complexity(m_spec::RidgeRegressionMachineSpec, m)
    params = fitted_params(m)
    c = 0.0
    if m_spec.include_bias
        c += abs(params.bias)^2
    end
    for (feature, coef) in params.coefs
        c += abs(coef)^2
    end
    return m_spec.lambda_ridge * c
end
