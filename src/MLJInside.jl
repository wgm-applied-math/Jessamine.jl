# For using MLJ machinery applied to the output of a genome

export LinearModelMachineSpec

using MLJLinearModels

@kwdef struct LinearModelMachineSpec <: AbstractMachineSpec
    model::Model
    lambda_ridge::Real
    lambda_parameter::Real
    lambda_operand::Real
    train_rows::AbstractVector
    test_rows::AbstractVector
end

"""
    machine_init(machine_spec::LinearModelMachineSpec, X, y)

Return `MLJ.machine(machine_spec.model, X, y)`
"""
function machine_init(m_spec::LinearModelMachineSpec, X, y)
    return MLJ.machine(m_spec.model, X, y)
end

"""
    machine_complexity(machine_spec::LinearModelMachineSpec, machine)

Return the sum of squares of the numbers in
`MLJ.fitted_params(m).coefs`
"""

function machine_complexity(m_spec::LinearModelMachineSpec, m)
    params = MLJ.fitted_params(m)
    c = 0.0

    # Hastie et al exclute the bias from the regularization term.
    # That way, addition of a constant to every given y
    # results in that same constant being added to the prediction.
    #     c += abs(params.bias)^2

    for (feature, coef) in params.coefs
        c += abs(coef)^2
    end
    return m_spec.lambda_ridge * c
end
