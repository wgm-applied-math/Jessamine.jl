# Make language server happy
export AbstractLinearModelResult
export BasicLinearModelResult, coefficients, intercept

export extend_if_singleton

"""
Abstract base type for the results of fitting linear models.
"""
abstract type AbstractLinearModelResult <: AbstractModelResult end

"""
A `BasicLinearModelResult` has a vector of coefficients
plust a bias or intercept.

"""
@kwdef struct BasicLinearModelResult <: AbstractLinearModelResult
    coefficients::Vector{Float64}
    intercept::Float64
end

BasicLinearModelResult(b) = BasicLinearModelResult(b, 0.0)

function coefficients(r::BasicLinearModelResult)
    return r.coefficients
end

function intercept(r::BasicLinearModelResult)
    return r.intercept
end

"""
    model_predict(mr::AbstractLinearModelResult, xs::AbstractArray{<:AbstractArray}; kw_args...)

Stack the columns `xs` into a matrix `X` and call `model_predict`.
The `kw_args` are splatted in.
"""
function model_predict(
        mr::AbstractLinearModelResult, xs::AbstractVector{<:AbstractVector}; kw_args...)
    X = stack(xs)
    return model_predict(mr, X; kw_args...)
end

function model_predict(
        mr::AbstractLinearModelResult, xs::NamedTuple; kw_args...)
    return model_predict(mr, fieldvalues(xs); kw_args...)
end

function model_predict(
        mr::AbstractLinearModelResult, xs::Tuple; kw_args...)
    return model_predict(mr, collect(xs); kw_args...)
end


"""
    model_predict(lmr::AbstractLinearModelResult, X::AbstractMatrix)

Return `X * coefficients(lmr) + intercept(lmr)`
"""
function model_predict(lmr::AbstractLinearModelResult, X::AbstractMatrix)
    return X * coefficients(lmr) .+ intercept(lmr)
end

"""
    model_predict(lmr::AbstractLinearModelResult, x::AbstractArray; kw_args...)

Return `dot(x, coefficients(lmr)) + intercept(lmr)`.
"""

function model_predict(
        lmr::AbstractLinearModelResult, x::AbstractArray;
        kw_args...
)
    # This is a bit of a hack to get the right type.
    # We need to convert the array to a matrix.
    return dot(x, coefficients(lmr)) + intercept(lmr)
end

"""
    extend_if_singleton(v::AbstractArray, m::Int)

If `v` is a singleton, as in `v = [x]`,
return `[x, x, ...]`, that is, a vector of `m` copies of `x`.
Otherwise, assert that `length(v) == m` and return `v`.
So the result is always a vector of length `m`.
"""
function extend_if_singleton(v::AbstractArray, m::Int)
    if length(v) == 1
        return fill(v[1], m)
    else
        @assert length(v)==m "m = $m, v = $v"
        return v
    end
end

function extend_if_singleton(v::AbstractArray, shape::Tuple{Int})
    (m,) = shape
    return extend_if_singleton(v, m)
end
