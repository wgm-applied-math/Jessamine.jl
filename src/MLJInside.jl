# For using MLJ machinery applied to the output of a genome

export LinearModelMachineSpec
export sum_sq_coefs

using MLJLinearModels

@kwdef struct LinearModelMachineSpec <: AbstractMachineSpec
    model::Model
    lambda_model::Float64
    lambda_parameter::Float64
    lambda_operand::Float64
    train_rows::AbstractVector
    performance::Any # callable
end
