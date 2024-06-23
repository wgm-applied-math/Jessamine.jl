export run_genome_symbolic

"""
    run_genome_symbolic(g_spec, genome; paramter_sym=:p, input_sym=:x)

Build a symbolic form for the output of the final time step of
running `genome`.  The parameter vector and input vector are
`Symbolics` objects of the form `p[j]` and `x[j]`.  The variable
names can be specified with the keyword arguments.

Return a named tuple `(p, x, w)` where `p` and `x`,
are vectors of `Symbolics` objects used to represent
genome parameters and inputs; and `w` is
a vector of genome outputs in symbolic form.
"""
function run_genome_symbolic(
        g_spec::GenomeSpec,
        genome::Genome;
        parameter_sym = :p,
        input_sym = :x)
    p = Symbolics.variables(parameter_sym, 1:(g_spec.parameter_size))
    x = Symbolics.variables(input_sym, 1:(g_spec.input_size))

    w = run_genome(g_spec, genome, p, x)[end][1:(g_spec.output_size)]
    return (p = p, x = x, w = w)
end
