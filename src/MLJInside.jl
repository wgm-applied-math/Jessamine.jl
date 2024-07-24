# For using MLJ machinery applied to the output of a genome
export BasicModelMachineSpec
export LinearModelMachineSpec
export sum_sq_coefs

using MLJLinearModels

@kwdef struct BasicModelMachineSpec <: AbstractMachineSpec
    model::Model
    lambda_parameter::Float64
    lambda_operand::Float64
    train_rows::AbstractVector
    performance::Any # callable
end

@kwdef struct LinearModelMachineSpec <: AbstractMachineSpec
    model::Model
    lambda_model::Float64
    lambda_parameter::Float64
    lambda_operand::Float64
    train_rows::AbstractVector
    performance::Any # callable
end

"""
    machine_complexity(mn_spec, m)

Numerical complexity of a machine.

Return the sum of squares of `MLJ.fitted_params(m)`
multiplied by `mn_spec.lambda_model`.
"""
function machine_complexity(mn_spec::LinearModelMachineSpec, m)
    params = MLJ.fitted_params(m)

    # Hastie et al exclude the bias from the regularization term.
    # That way, addition of a constant to every given y
    # results in that same constant being added to the prediction.
    #     c += abs(params.intercept)^2
    return mn_spec.lambda_model * sum_sq_coefs(params.coefs)
end

function sum_sq_coefs(coefs::Dict)
    c = 0.0
    for (feature, coef) in coefs
        c += abs(coef)^2
    end
    return c
end

function sum_sq_coefs(coefs::AbstractVector{<:Pair})
    c = 0.0
    for (feature, coef) in coefs
        c += abs(coef)^2
    end
    return c
end

function sum_sq_coefs(coefs::AbstractVector{<:Number})
    return m_spec.lambda_model * dot(coefs, coefs)
end
