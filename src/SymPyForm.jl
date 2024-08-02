# Make language server happy
if false
    using SymPy
    include("GenomeCore.jl")
end

export eval_time_step_sympy, show_sympy, run_genome_sympy
export replace_near_integer

"""
    eval_time_step_sympy(g_spec, genome; output_sym, scratch_sym, parameter_sym, input_sym)

Return symbolic objects representing a single time step of `genome`.
The result is a named tuple with fields
`z`, `t`, `p`, `x` for the `SymPy` objects for those variables;
`c` and `c_next` for the current and future cell state in symbolic form.
"""

function eval_time_step_sympy(
        g_spec::GenomeSpec,
        genome::Genome;
        output_sym = :z,
        scratch_sym = :t,
        parameter_sym = :p,
        input_sym = :x)
    z = [Sym("$output_sym$j") for j in 1:(g_spec.output_size)]
    t = [Sym("$scratch_sym$j") for j in 1:(g_spec.scratch_size)]
    p = [Sym("$parameter_sym$j") for j in 1:(g_spec.parameter_size)]
    x = [Sym("$input_sym$j") for j in 1:(g_spec.input_size)]
    c = CellState(z, t, p, x)
    c_next = eval_time_step(c, genome)
    return (z = z, t = t, p = p, x = x, c = c, c_next = c_next)
end

function show_sympy(
        g_spec::GenomeSpec,
        genome::Genome;
        output_sym = :z,
        scratch_sym = :t,
        parameter_sym = :p,
        input_sym = :x)
    sym_res = eval_time_step_sympy(
        g_spec,
        genome;
        output_sym = output_sym,
        scratch_sym = scratch_sym,
        parameter_sym = parameter_sym,
        input_sym = input_sym)
    return hcat(
        1:length(sym_res.c),
        flat_workspace(sym_res.c),
        flat_workspace(sym_res.c_next))
end

function run_genome_sympy(
        g_spec::GenomeSpec,
        genome::AbstractGenome;
        parameter_sym = :p,
        input_sym = :x)
    p = [Sym("$parameter_sym$j") for j in 1:(g_spec.parameter_size)]
    x = [Sym("$input_sym$j") for j in 1:(g_spec.input_size)]
    z = run_genome(g_spec, genome, p, x)[end]
    return (p = p, x = x, z = z)
end

function replace_near_integer(expr::Sym; tolerance = 1.0e-10)
    if expr.func.𝑓.__name__ == "Float"
        x = expr.evalf()
        k = round(x)
        if k < typemax(Int) && abs(expr - k) < tolerance
            return Sym(k)
        else
            return expr
        end
    else
        new_args = map(expr.args) do expr
            replace_near_integer(expr; tolerance = tolerance)
        end
        return expr.func(new_args)
    end
end
