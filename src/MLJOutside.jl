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
    include("SymbolicForm.jl")
end

export EpochSpec

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

@kwdef mutable struct EpochSpec
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

# I really need something like
# struct JModel{Super} <: Super   inner_model::Super ... end
# so that JModel{Deterministic} <: Deterministic.
# However, this is not possible in Julia's type system.
# The syntax struct Thing{T} <: T ... end
# results in Thing{T} <: (T where T),
# or equivalently Thing{T} <: Any.
# Thing{T} <: T evaluates to false.

# So we resort to macros.

# Complication: @mlj_model only works on a
# non-parametric struct.

macro declareJessamineModel(model_type, default_inner_model)
    struct_name = Symbol("Jessamine$(model_type)")
    return quote
        export $(esc(struct_name))
        @mlj_model mutable struct $(esc(struct_name)) <: $model_type

            inner_model::$model_type = $default_inner_model

            rng::AbstractRNG = Random.default_rng()
            init_arity_dist::Distribution = DiscreteNonParametric(
                [1, 2, 3], [0.25, 0.5, 0.25])

            # *MachineModelSpec
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

            neighborhoods::AbstractVector{EpochSpec} =
                default_neighborhoods::(length(_) >= 1)
            num_epochs::Int = 10::(_ >= 0)
            simplifier::Union{Nothing, EpochSpec} = default_simplifier
            sense::Any = MinSense
            stop_threshold::Union{Nothing, <:Number} = 0.001

            logging_generation_modulus::Int = 10::(_ >= 1)
        end
    end
end

@declareJessamineModel Deterministic RidgeRegressor(lambda = 0.01)
@declareJessamineModel Probabilistic LogisticClassifier()

"""
    machine_spec_type(t_MLJ_model::Type)::Type{<:AbstractMachineSpec}

Return the Jessamine type (subtype of `AbstractMachineSpec`)
corresponding to the given MLJ model type.
The default is to return `BasicModelMachineSpec`
"""
function machine_spec_type(::Type)
    return BasicModelMachineSpec
end

"""
    machine_spec_type(::Type{RidgeRegressor})

Return `LinearModelMachineSpec`
"""
function machine_spec_type(::Type{RidgeRegressor})
    return LinearModelMachineSpec
end

"""
    machine_spec_type(::Type{LogisticClassifier})

Return `LinearModelMachineSpec`
"""
function machine_spec_type(::Type{LogisticClassifier})
    return LinearModelMachineSpec
end



"""
A union of `JessamineDeterministic` and `JessamineProbabilistic`
"""

const JessamineModel = Union{JessamineDeterministic, JessamineProbabilistic}

function build_specs!(
    jm::JessamineModel,
    X,
    y,
    verbosity)
    xs = Tables.columns(X)
    jm.input_size = length(xs)
    g_spec = convert(GenomeSpec, jm)
    mn_spec = convert(machine_spec_type(jm.model), jm)
    function grow_and_rate(rng, g_spec, genome)
        return machine_grow_and_rate(
            xs, y, g_spec, genome, mn_spec,
            DefaultSolverSpec())
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
        convert(MutationSpec, jm.simplifier),
        convert(SelectionSpec, jm.simplifier),
        grow_and_rate,
        jm.simplifier.num_generations
    )

    return (g_spec=g_spec,
            mn_spec=mn_spec,
            neighborhoods=neighborhoods,
            simp_spec=simp_spec)
end

function MLJModelInterface.fit(
    jm::JessamineModel,
    verbosity::Int,
    X,
    y)
    specs = build_specs!(jm, X, y, verbosity)
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
        generation_mod = jm.logging_generation_modulus,
        verbosity = verbosity)
    pop_simp = evolution_loop(
        jm.rng,
        specs.simp_spec,
        pop_next;
        sense = jm.sense,
        stop_threshold = jm.stop_threshold,
        generation_mod = jm.logging_generation_modulus,
        verbosity = verbosity)
    best_agent = pop_simp.agents[1]
    symbolic_result = run_genome_symbolic(specs.g_spec, best_agent.genome)
    fit_result = (
        g_spec = specs.g_spec,
        best_agent = best_agent,
        inner_fitted_params = MLJ.fitted_params(best_agent.extra.m))
    report = (
        rating = best_agent.rating,
        symbolic_result = symbolic_result)
    # This way you can continue with update() I think...
    cache = pop_simp
    return (fit_result, report, cache)
end
