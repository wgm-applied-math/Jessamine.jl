export splat_or_default
export Add, Multiply

"""
    splat_or_default(op, def)

Return a function for which `f([]) = def` and `f([x1...]) = op(x1...)`.
The `Base.splat` function doesn't have a way to deal with the first case.
"""
function splat_or_default(op, def)
    function f(operands)
        if isempty(operands)
            return def
        else
            return splat(op)(operands)
        end
    end
    return f
end

"Add elements of the work space vector."
struct Add <: AbstractGeneOp end
as_function(::Add) = splat_or_default(.+, 0)
short_show(io::IO, ::Add) = print(io, "add")

"Multiply elements of the work space vector."
struct Multiply <: AbstractGeneOp end
as_function(::Multiply) = splat_or_default(.*, 1)
short_show(io::IO, ::Multiply) = print(io, "mul")
