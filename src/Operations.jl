# Make language server happy
if false
    include("GenomeCore.jl")
end

abstract type AbstractUnaryOp end

abstract type AbstractMultiOp <: AbstractGeneOp end

export splat_or_default
export AbstractUnaryOp, AbstractMultiOp
export Add, Multiply, Subtract
export UnaryComposition
export ReciprocalMultiply, ReciprocalAdd, ReciprocalSubtract
export FzAnd, FzOr, FzNand, FzNor
export Maximum, Minimum
export Sign, SignAdd, SignSubtract
export Sigmoid, SigmoidAdd, SigmoidSubtract

"""
    splat_or_default(op, def, operands)

Return the result of applying `op` to the operands,
with `op([]) = def` and `op([x1...]) = op(x1...)`.
The `Base.splat` function doesn't have a way to deal with the first case.
"""
function splat_or_default(op, def, operands)
    if isempty(operands)
        return def
    else
        return splat(op)(operands)
    end
end

"Add operands."
struct Add <: AbstractMultiOp end

short_show(io::IO, ::Add) = print(io, "add")

op_eval(::Add, operands) = splat_or_default(.+, 0.0, operands)

function op_eval_add_into!(dest::AbstractVector, ::Add, operands::AbstractVector)
    for operand in operands
        dest .= dest .+ operand
    end
end

function to_expr(::Add, cs, operands)
    if isempty(operands)
        return [0.0]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        return Expr(:call, :.+, (:($cs.$field[$j]) for (field, j) in operands)...)
    end
end

"Subtract LISP-style: `0 + x[1] - x[2] - x[3]...`"
struct Subtract <: AbstractMultiOp end

short_show(io::IO, ::Subtract) = print(io, "sub")

function op_eval(::Subtract, operands)
    if isempty(operands)
        return 0.0
    else
        return operands[1] .- op_eval(Add(), operands[2:end])
    end
end

function op_eval_add_into!(dest::AbstractVector, ::Subtract, operands::AbstractVector)
    if !isempty(operands)
        dest .= dest .+ operands[1]
        for operand in operands[2:end]
            dest .= dest .- operand
        end
    end
end

function to_expr(::Subtract, cs, operands)
    if isempty(operands)
        return [0.0]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        field1, j1 = operands[1]
        return Expr(:call, :.-, :($cs.$field1[$j1]),
            Expr(:call, :.+, (:($cs.$field[$j]) for (field, j) in operands[2:end])...))
    end
end

"Multiply operands."
struct Multiply <: AbstractMultiOp end

short_show(io::IO, ::Multiply) = print(io, "mul")

op_eval(::Multiply, operands) = splat_or_default(.*, 1.0, operands)

function to_expr(::Multiply, cs, operands)
    if isempty(operands)
        return [1.0]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        return Expr(:call, :.*, (:($cs.$field[$j]) for (field, j) in operands)...)
    end
end

"Do a multi-arity operation and apply a unary operation"
@kwdef struct UnaryComposition{Un <: AbstractUnaryOp, Multi <: AbstractMultiOp} <:
              AbstractMultiOp
    unary::Un = Un()
    multi::Multi = Multi()
end

function short_show(io::IO, c::UnaryComposition)
    short_show(io, c.unary)
    print(io, "@")
    short_show(io, c.multi)
end

function op_eval(c::UnaryComposition, operands)
    return un_op_eval(c.unary, op_eval(c.multi, operands))
end

function to_expr(c::UnaryComposition, cs, operands)
    inner_expr = to_expr(c.multi, cs, operands)
    outer_expr = to_expr(c.unary, inner_expr)
    return outer_expr
end

struct Reciprocal <: AbstractUnaryOp end
short_show(io::IO, ::Reciprocal) = print(io, "rcp")
un_op_eval(::Reciprocal, t) = 1.0 ./ t
to_expr(::Reciprocal, expr) = :(1.0 ./ $expr)

"Multiply operands and return the reciprocal."
const ReciprocalMultiply = UnaryComposition{Reciprocal, Multiply}

"Add operands and return the reciprocal."
const ReciprocalAdd = UnaryComposition{Reciprocal, Add}

"Subtract operands and return the reciprocal."
const ReciprocalSubtract = UnaryComposition{Reciprocal, Subtract}

"Return fuzzy AND of the operands"
struct FzAnd <: AbstractMultiOp end

short_show(io::IO, ::FzAnd) = print(io, "fzAnd")

op_eval(::FzAnd, operands) = splat_or_default(.*, 1.0, operands)

function to_expr(::FzAnd, cs, operands)
    if isempty(operands)
        return [1.0]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        return Expr(:call, :.*, (:($cs.$field[$j]) for (field, j) in operands)...)
    end
end

"Return fuzzy OR of the operands."
struct FzOr <: AbstractMultiOp end

short_show(io::IO, ::FzOr) = print(io, "fzOr")

function op_eval(::FzOr, operands)
    return 1.0 .- op_eval(FzNor(), operands)
end

function to_expr(::FzOr, cs, operands)
    return Expr(:call, :.-, 1.0, to_expr(FzNor(), cs, operands))
end

"Return fuzzy NAND of the operands."
struct FzNand <: AbstractMultiOp end

short_show(io::IO, ::FzNand) = print(io, "fzNand")

function op_eval(::FzNand, operands)
    return 1.0 .- op_eval(FzAnd(), operands)
end

function to_expr(::FzNand, cs, operands)
    return Expr(:call, :.-, 1.0, to_expr(FzAnd(), cs, operands))
end

"Return fuzzy NOR of the operands."
struct FzNor <: AbstractMultiOp end

short_show(io::IO, ::FzNor) = print(io, "fzNor")

function op_eval(::FzNor, operands)
    return splat_or_default(.*, 1.0, 1.0 .- operand for operand in operands)
end

function to_expr(::FzNor, cs, operands)
    if isempty(operands)
        return [1.0]
    elseif length(operands) == 1
        field, j = operands[1]
        return :(1.0 .- $cs.$field[$j])
    else
        return Expr(:call, :.*, (:(1.0 .- $cs.$field[$j]) for (field, j) in operands)...)
    end
end

"Return the sign of the operand"
struct Sign <: AbstractUnaryOp end
short_show(io::IO, ::Sign) = print(io, "sign")
un_op_eval(::Sign, t) = sign.(t)
to_expr(::Sign, expr) = :(sign.($expr))

"Return the sign of the sum of the operands"
const SignAdd = UnaryComposition{Sign, Add}

"Return the sign of the difference of the operands"
const SignSubtract = UnaryComposition{Sign, Subtract}

"Apply the exponential sigmoid to the operand"
struct Sigmoid <: AbstractUnaryOp end
short_show(io::IO, ::Sign) = print(io, "sigmoid")
un_op_eval(::Sigmoid, t) = 1 ./ (1 .+ exp.(-t))
to_expr(::Sigmoid, expr) = :(1 ./ (1 .+ exp.(-$expr)))

"Apply the exponential sigmoid to the sum of the operands"
const SigmoidAdd = UnaryComposition{Sigmoid, Add}

"Apply the exponential sigmoid to the difference of the operands"
const SigmoidSubtract = UnaryComposition{Sigmoid, Subtract}

"Return the maximum of the operands."
struct Maximum <: AbstractMultiOp end

short_show(io::IO, ::Maximum) = print(io, "max")

function op_eval(::Maximum, operands)
    if isempty(operands)
        return -Inf
    elseif length(operands) == 1
        return operands[1]
    else
        return max.(operands...)
    end
end

function to_expr(::Maximum, cs, operands)
    if isempty(operands)
        return [-Inf]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        return quote
            max.($((:($cs.$field[$j]) for (field, j) in operands)...))
        end
    end
end

"Return the minimum of the operands."
struct Minimum <: AbstractMultiOp end

short_show(io::IO, ::Minimum) = print(io, "min")

function op_eval(::Minimum, operands)
    if isempty(operands)
        return -Inf
    elseif length(operands) == 1
        return operands[1]
    else
        return min.(operands...)
    end
end

function to_expr(::Minimum, cs, operands)
    if isempty(operands)
        return [Inf]
    elseif length(operands) == 1
        field, j = operands[1]
        return :($cs.$field[$j])
    else
        return quote
            min.($((:($cs.$field[$j]) for (field, j) in operands)...))
        end
    end
end
