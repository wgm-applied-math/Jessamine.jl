export splat_or_default
export Add, Multiply, ReciprocalMultiply

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

"Add elements of the workspace vector."
struct Add <: AbstractGeneOp end
op_eval(::Add, operands) = splat_or_default(.+, 0, operands)
short_show(io::IO, ::Add) = print(io, "add")

"Multiply elements of the workspace vector."
struct Multiply <: AbstractGeneOp end
op_eval(::Multiply, operands) = splat_or_default(.*, 1, operands)
short_show(io::IO, ::Multiply) = print(io, "mul")

"Multiply elements of the workspace vector and return the reciprocal."
struct ReciprocalMultiply <: AbstractGeneOp end
op_eval(::ReciprocalMultiply, operands) = 1 ./ splat_or_default(.*, 1, operands)
short_show(io::IO, ::ReciprocalMultiply) = print(io, "rcpm")