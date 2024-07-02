export splat_or_default
export Add, Multiply, ReciprocalMultiply
export Subtract

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
struct Add <: AbstractGeneOp end
short_show(io::IO, ::Add) = print(io, "add")
op_eval(::Add, operands) = splat_or_default(.+, 0, operands)

"Subtract LISP-style: `0 + x[1] - x[2] - x[3]...`"
struct Subtract <: AbstractGeneOp end
short_show(io::IO, ::Subtract) = print(io, "sub")
function op_eval(::Subtract, operands)
    if isempty(operands)
        return 0
    else
        return operands[1] - sum(operands[2:end])
    end
end

"Multiply operands."
struct Multiply <: AbstractGeneOp end
short_show(io::IO, ::Multiply) = print(io, "mul")
op_eval(::Multiply, operands) = splat_or_default(.*, 1, operands)


"Multiply operands and return the reciprocal."
struct ReciprocalMultiply <: AbstractGeneOp end
short_show(io::IO, ::ReciprocalMultiply) = print(io, "rcpm")
op_eval(::ReciprocalMultiply, operands) = 1 ./ splat_or_default(.*, 1, operands)
