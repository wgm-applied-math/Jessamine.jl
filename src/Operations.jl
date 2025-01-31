# !!! TODO Look into NaNMath.jl which returns NaN instead of throwing an error for things
# like sqrt(negative)

# Make language server happy
if false
    # include("GenomeCore.jl")
end

abstract type AbstractUnaryOp end

abstract type AbstractMultiOp <: AbstractGeneOp end

export splat_or_default
export AbstractUnaryOp, AbstractMultiOp
export Add, Multiply, Subtract
export UnaryComposition
export Reciprocal, ReciprocalMultiply, ReciprocalAdd, ReciprocalSubtract
export FzAnd, FzOr, FzNand, FzNor
export Maximum, Minimum
export Sign, SignAdd, SignSubtract
export Sigmoid, SigmoidAdd, SigmoidSubtract
export define_unary_op

"""
    splat_or_default(op, def, workspace, indices)

Return the result of applying `op` to the `operands = workspace[indices]`,
with `op([]) = def` and `op([x1...]) = op(x1...)`.
The `Base.splat` function doesn't have a way to deal with the first case.
"""
function splat_or_default(op, def, operands)
    if isempty(operands)
        return def
    else
        # Splatting is oddly slow.
        # Specialized implementations are given for more specific types.
        return op(operands...)
    end
end

function splat_or_default(op, def, workspace::AbstractVector{<:Number}, indices::AbstractVector)
    if isempty(indices)
        return def
    else
        res = workspace[indices[1]]
        for r in indices[2:end]
            res = op(res, workspace[r])
        end
        return res
    end
end


function splat_or_default(op, def, workspace::AbstractVector{<:AbstractVector}, indices::AbstractVector)
    if isempty(indices)
        return def
    else
        n = 0
        for v in workspace
            n = max(n, length(v))
        end
        e_type = eltype(workspace[1])
        res = Vector{e_type}(undef, n)
        res .= workspace[indices[1]]
        for j in 2:length(indices)
            res .= op(res, workspace[indices[j]])
        end
        return res
    end
end


"Add operands."
struct Add <: AbstractMultiOp end

short_show(io::IO, ::Add) = print(io, "add")

op_eval(::Add, workspace, indices) = splat_or_default(.+, 0.0, workspace, indices)

function op_eval_add_into!(dest::AbstractArray, ::Add, workspace::AbstractArray, indices::AbstractArray)
    for j in indices
        dest .+= workspace[j]
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

function op_eval(::Subtract, workspace, indices)
    if isempty(indices)
        return 0.0
    else
        return workspace[indices[1]] .- op_eval(Add(), workspace, indices[2:end])
    end
end

function op_eval_add_into!(dest::AbstractArray, ::Subtract, workspace::AbstractArray, indices::AbstractArray)
    if !isempty(indices)
        dest .+= workspace[indices[1]]
        for j in indices[2:end]
            dest .-= workspace[j]
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

op_eval(::Multiply, workspace, indices) = splat_or_default(.*, 1.0, workspace, indices)

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

function op_eval(c::UnaryComposition, workspace, indices)
    return un_op_eval(c.unary, op_eval(c.multi, workspace, indices))
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

op_eval(::FzAnd, workspace, indices) = splat_or_default(.*, 1.0, workspace, indices)

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

function op_eval(::FzOr, workspace, indices)
    return 1.0 .- op_eval(FzNor(), workspace, indices)
end

function to_expr(::FzOr, cs, operands)
    return Expr(:call, :.-, 1.0, to_expr(FzNor(), cs, operands))
end

"Return fuzzy NAND of the operands."
struct FzNand <: AbstractMultiOp end

short_show(io::IO, ::FzNand) = print(io, "fzNand")

function op_eval(::FzNand, workspace, indices)
    return 1.0 .- op_eval(FzAnd(), workspace, indices)
end

function to_expr(::FzNand, cs, operands)
    return Expr(:call, :.-, 1.0, to_expr(FzAnd(), cs, operands))
end

"Return fuzzy NOR of the operands."
struct FzNor <: AbstractMultiOp end

short_show(io::IO, ::FzNor) = print(io, "fzNor")

function op_eval(::FzNor, workspace, indices)
    return splat_or_default((acc, x) -> acc .* (1.0 .- x), 1.0, workspace, indices)
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
short_show(io::IO, ::Sigmoid) = print(io, "sigmoid")
un_op_eval(::Sigmoid, t) = 1 ./ (1 .+ exp.(-t))
to_expr(::Sigmoid, expr) = :(1 ./ (1 .+ exp.(-$expr)))

"Apply the exponential sigmoid to the sum of the operands"
const SigmoidAdd = UnaryComposition{Sigmoid, Add}

"Apply the exponential sigmoid to the difference of the operands"
const SigmoidSubtract = UnaryComposition{Sigmoid, Subtract}

"Return the maximum of the operands."
struct Maximum <: AbstractMultiOp end

short_show(io::IO, ::Maximum) = print(io, "max")

function op_eval(::Maximum, workspace, indices)
    if isempty(indices)
        return -Inf
    elseif length(indices) == 1
        return workspace[indices[1]]
    else
        # This can be sped up...
        return max.(workspace[indices]...)
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

function op_eval(::Minimum, workspace, indices)
    if isempty(indices)
        return -Inf
    elseif length(indices) == 1
        return workspace[indices[1]]
    else
        # This can be sped up...
        return min.(workspace[indices]...)
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

macro define_unary_op(struct_name, function_name)
    return quote
        export $struct_name
        export $(Symbol(string(struct_name) * "Add"))
        export $(Symbol(string(struct_name) * "Subtract"))
        export $(Symbol(string(struct_name) * "Multiply"))
        @doc "Return "*string($function_name)*" of the operand"
        struct $struct_name <: AbstractUnaryOp end
        global short_show
        short_show(io::IO, ::($struct_name)) = print(io, string($function_name))
        global un_op_eval
        un_op_eval(::($struct_name), t) = ($function_name).(t)
        global to_expr
        to_expr(::($struct_name), expr) = :($($function_name).($expr))
        @doc "Return "*string($function_name)*" applied to the sum of the operands"
        const $(Symbol(string(struct_name) * "Add")) = UnaryComposition{$struct_name, Add}
        @doc "Return "*string($function_name)*" applied to the result of subtraction of the operands"
        const $(Symbol(string(struct_name) * "Subtract")) = UnaryComposition{
            $struct_name, Subtract}
        @doc "Return "*string($function_name)*" applied to the product of the operands"
        const $(Symbol(string(struct_name) * "Multiply")) = UnaryComposition{
            $struct_name, Multiply}
    end
end

@define_unary_op Sqrt sqrt
@define_unary_op Exp exp
@define_unary_op Log log
@define_unary_op Sin sin
@define_unary_op Cos cos
@define_unary_op Tan tan
@define_unary_op ASin asin
@define_unary_op ACos acos
@define_unary_op ATan atan
@define_unary_op Sinh sinh
@define_unary_op Cosh cosh
@define_unary_op Tanh tanh
@define_unary_op ASinh asinh
@define_unary_op ACosh acosh
@define_unary_op ATanh atanh

export PolynomialInventory
const PolynomialInventory = [Add(), Subtract(), Multiply()]

export RationalFunctionInventory
const RationalFunctionInventory = vcat(
    PolynomialInventory, [ReciprocalAdd(), ReciprocalSubtract(), ReciprocalMultiply()])

export ExpLogInventory
const ExpLogInventory = vcat(RationalFunctionInventory,
    reshape(
        [UnaryComposition{un_op, bin_op}()
         for un_op in [Sqrt, Exp, Log], bin_op in [Add, Subtract, Multiply]],
        :))

export TrigInventory
const TrigInventory = vcat(ExpLogInventory,
    reshape(
        [UnaryComposition{un_op, bin_op}()
         for un_op in [Sin, Cos, Tan, ASin, ACos, ATan], bin_op in [Add, Subtract, Multiply]],
        :))

export HyperbolicInventory
const HyperbolicInventory = vcat(ExpLogInventory,
        reshape(
            [UnaryComposition{un_op, bin_op}()
             for un_op in [Sinh, Cosh, Tanh, ASinh, ACosh, ATanh], bin_op in [Add, Subtract, Multiply]],
            :))
