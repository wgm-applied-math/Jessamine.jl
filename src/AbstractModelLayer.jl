export AbstractModelResult
export model_predict
export model_basic_symbolic_output, model_symbolic_output
export model_basic_sympy_output, model_sympy_output

"""
Abstract base type for the results of fitting a model
to data.  This object can then be used for predictions.
"""
abstract type AbstractModelResult end

"""
    model_predict(mr::AbstractModelResult, X; kw_args...)

Given the table `X` where the columns are inputs (predictors)
and the rows are points, use `mr` to predict the target value
for all the points.
"""
function model_predict end

"""
    model_predict(g_spec::GenomeSpec, agent::Agent, xs; kw_args...)

Run `agent.genome` on inputs `xs` and `agent.parameter`, and
form the linear combination of the genome's
outputs using the coefficients `agent.extra`.
"""
function model_predict(
        g_spec::GenomeSpec,
        agent::Agent{<:Number, <:AbstractGenome, <:AbstractArray, <:AbstractModelResult},
        x_table;
        kw_args...
)
    xs = separate_columns(x_table)
    num_rows = length(xs[1])
    last_round = run_genome_to_last(g_spec, agent.genome, agent.parameter, xs)
    data_cols = map(u -> extend_if_singleton(u, num_rows), last_round)
    output_col_names = map(1:(g_spec.output_size)) do t
        "z$t"
    end
    Z_table = namedtuple(output_col_names, data_cols)
    return model_predict(agent.extra, Z_table; kw_args...)
end
