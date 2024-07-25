# Make language server happy
if false
    using MLJ
    using MLJModelInterface
    using Optimization
    using Random
    using Tables
    include("Evolution.jl")
    include("GenomeCore.jl")
    include("Operations.jl")
    include("MachineLayer.jl")
    include("MLJInside.jl")
    include("Mutation.jl")
end

export EpochSpec
export JessamineDeterministicRegressor

# This is kind of a mess because in MLJ, the input size
# is not really available until the call to machine().
# Which means that the GenomeSpec can't be part of the model.
# Which means that EvolutionSpec can't be part of the model.
# Hence the following massive config system.
#
# This is why there are overloads for convert() for spec types.
# convert(::Type{SpecType}, x) constructs the spec object using
# fields of x with the same name as a field in the SpecType.
# That way we can call convert(MutationSpec, my_epoc_spec) for
# example.

@kwdef struct EpochSpec
    # MutationSpec
    p_mutate_op::Float64 = 0.1
    p_mutate_index::Float64 = 0.1
    p_duplicate_index::Float64 = 0.01
    p_delete_index::Float64 = 0.01
    p_duplicate_instruction::Float64 = 0.001
    p_delete_instruction::Float64 = 0.001
    op_inventory::Vector{AbstractGeneOp} = [Add(), Subtract(), Multiply()]
    op_probabilities::Union{Vector{Float64}, Nothing} = [0.25, 0.25, 0.5]

    # SelectionSpec
    num_to_keep::Int = 20
    num_to_generate::Int = 40
    p_take_better::Float64 = 0.8
    p_take_very_best::Float64 = 0.25

    # EvolutionSpec
    num_generations::Int = 10
    stop_on_innovation::Bool = false
end


# Set up default neighborhood structure

e_spec_1 = EpochSpec()
e_spec_2 = EpochSpec(
    p_mutate_op=0.2,
    p_mutate_index=0.2,
    p_duplicate_index=0.02,
    p_delete_index=0.02,
    p_duplicate_instruction=0.002,
    p_delete_instruction=0.002,
    num_generations=40,
    stop_on_innovation=true)
e_spec_3 = EpochSpec(
    p_mutate_op=0.2,
    p_mutate_index=0.2,
    p_duplicate_index=0.02,
    p_delete_index=0.02,
    p_duplicate_instruction=0.002,
    p_delete_instruction=0.002,
    num_to_keep=10,
    num_to_generate=50,
    p_take_better=0.6,
    p_take_very_best=0.0,
    num_generations=40,
    stop_on_innovation=true)
e_spec_4 = EpochSpec(
    p_mutate_op=0.3,
    p_mutate_index=0.3,
    p_duplicate_index=0.03,
    p_delete_index=0.03,
    p_duplicate_instruction=0.003,
    p_delete_instruction=0.003,
    num_to_keep=10,
    num_to_generate=50,
    p_take_better=0.6,
    p_take_very_best=0.0,
    num_generations=40,
    stop_on_innovation=true)
e_spec_5 = EpochSpec(
    p_mutate_op=0.3,
    p_mutate_index=0.3,
    p_duplicate_index=0.03,
    p_delete_index=0.03,
    p_duplicate_instruction=0.003,
    p_delete_instruction=0.003,
    num_to_keep=10,
    num_to_generate=50,
    p_take_better=0.5,
    p_take_very_best=0.0,
    num_generations=40,
    stop_on_innovation=true)

const default_neighborhoods = [
    e_spec_1, e_spec_2, e_spec_3, e_spec_4, e_spec_5
]

const default_simplifier = EpochSpec(
    p_mutate_op=0.0,
    p_mutate_index=0.0,
    p_duplicate_index=0.0,
    p_delete_index=0.1,
    p_duplicate_instruction=0.0,
    p_delete_instruction=0.1,
    num_generations=100
)

const default_deterministic_regressor_inner_model = RidgeRegressor(lambda = 0.01)

@mlj_model mutable struct JessamineDeterministicRegressor{Inner<:Deterministic}
    <: Deterministic
    inner_model::Inner = default_deterministic_regressor_inner_model

    rng::AbstractRNG = Random.default_rng()
    init_arity_dist::Distribution = DiscreteNonParametric([1, 2, 3], [0.25, 0.5, 0.25])

    # *MachineModelSpec
    ModelMachineSpec::Type{<:AbstractMachineSpec} = LinearModelMachineSpec
    lambda_model::Float64 = 0.01::(_ >= 0.0)
    lambda_parameter::Float64 = 0.01::(_ >= 0.0)
    lambda_operand::Float64 = 0.01::(_ >= 0)
    performance::Any = MLJ.l2

    # GenomeSpec
    output_size::Int = 6::(_ >= 0)
    scratch_size::Int = 6::(_ >= 0)
    parameter_size::Int = 2::(_ >= 0)
    input_size::Int = 0::(_ >= 0)
    num_time_steps::Int = 3::(_ >= 0)

    neighborhoods::AbstractVector{EpochSpec} = default_neighborhoods::(length(_) >= 1)
    num_epochs::Int = 10::(_ >= 0)
    simplifier::Union{Nothing,EpochSpec} = default_simplifier
    sense::Any = MinSense
    stop_threshold::Union{Nothing,<:Number} = 0.001

    log_generation_mod::Int = 10::(_ >= 1)
end

function build_specs!(
    jm::JessamineDeterministicRegressor,
    X,
    y,
    fit_kw_args)
    xs = Tables.columns(X)
    jm.input_size = length(xs)
    g_spec = convert(GenomeSpec, jm)
    mn_spec = convert(jm.ModelMachineSpec, jm)
    function grow_and_rate(rng, g_spec, genome)
        return machine_grow_and_rate(
            xs, y, g_spec, genome, mn_spec,
            DefaultSolverSpec(),
            fit_kw_args=fit_kw_args)
    end
    neighborhoods = map(jm.neighborhoods) do epoch_spec
        m_spec = convert(MutationSpec, epoch_spec)
        s_spec = convert(SelectionSpec, epoch_spec)
        e_spec = EvolutionSpec(
            g_spec,
            m_spec,
            s_spec,
            grow_and_rate,
            epoch_spec.num_generations,
            epoch_spec.stop_on_innovation)
        return e_spec
    end
    simp_spec = EvolutionSpec(
        g_spec,
        convert(MutationSpec, jm.simplifier)
        convert(SelectionSpec, jm.simplifier)
        grow_and_rate,
        jm.simplifier.num_generations
    )

    return (g_spec=g_spec,
            mn_spec=mn_spec,
            neighborhoods=neighborhoods,
            simp_spec=simp_spec)
end

function MLJModelInterface.fit(
    jm::JessamineDeterministicRegressor,
    verbosity::Int,
    X,
    y)
    specs = build_specs!(jm, X, y)
    pop_init = random_initial_population(
        jm.rng,
        specs.neighborhoods[1],
        jm.init_arity_dist,
        specs.grow_and_rate;
        sense = jm.sense)
    pop_next = vns_evolution_loop(
        jm.rng,
        specs.neighborhoods,
        jm.num_epochs,
        pop_init;
        stop_threshold = jm.stop_threshold,
        generation_mod = jm.log_generation_mod,
        enable_logging = (verbosity > 0))
    pop_simp = evolution_loop(
        jm.rng,
        specs.simp_spec,
        pop_next;
        sense = jm.sense,
        stop_threshold = jm.stop_threshold,
        generation_mod = jm.log_generation_mod,
        enable_logging = (verbosity > 0))
    best_agent = pop_simp.agents[1]
    sym_res = run_genome_symbolic(g_spec, best_agent.genome)
    fit_result = (
        g_spec = g_spec,
        best_agent = best_agent)
    report = (
        rating = best_agent.rating,
        sym_res = sym_res)
    cache = nothing
    return (fit_result, report, cache)
end
