export MutationSpec, MutationDist, RandomDuplicateDelete
export mutate

"""
A collection of parameters specifying the probabilities of various mutations.
"""
@kwdef struct MutationSpec
    "Probability that an operation changes"
    p_mutate_op::Float64 = 1.0e-3

    "Probability that an index changes"
    p_mutate_index::Float64 = 1.0e-3

    "Probability that an index is duplicated"
    p_duplicate_index::Float64 = 1.0e-4

    "Probability that an index is deleted"
    p_delete_index::Float64 = 1.0e-4

    "Probability that an instruction is duplicated"
    p_duplicate_instruction::Float64 = 1.0e-5

    "Probability that an instruction is deleted"
    p_delete_instruction::Float64 = 1.0e-5

    "Which operations are available"
    op_inventory::Vector{AbstractGeneOp} = []
end

"""
A collection of distribution objects built from a `MutationSpec` and used to produce random numbers determining mutations.
"""
@kwdef struct MutationDist
    d_mutate_op::Bernoulli
    d_mutate_index::Bernoulli
    d_new_src_index::DiscreteUniform
    d_duplicate_index::Bernoulli
    d_delete_index::Bernoulli
    d_duplicate_instruction::Bernoulli
    d_delete_instruction::Bernoulli
    op_inventory::Vector{AbstractGeneOp}
end

"""
    MutationDist(ms::MutationSpec, src_index_max::Int)

Use `ms` to build a `MutationDist`.
Specify that an operator index must be <= `src_index_max`.
"""
function MutationDist(ms::MutationSpec, src_index_max::Int)
    return MutationDist(
        Bernoulli(ms.p_mutate_op),
        Bernoulli(ms.p_mutate_index),
        DiscreteUniform(1, src_index_max),
        Bernoulli(ms.p_duplicate_index),
        Bernoulli(ms.p_delete_index),
        Bernoulli(ms.p_duplicate_instruction),
        Bernoulli(ms.p_delete_instruction),
        ms.op_inventory
    )
end

"""
    mutate(rng::AbstractRNG, md::MutationDist, x)

Randomly change `x`, using random numbers supplied by RNG and
the distributions specified by `md`.
"""
function mutate end

function mutate(rng::AbstractRNG, md::MutationDist, j::Int)
    if rand(rng, md.d_mutate_index)
        return rand(rng, md.d_new_src_index)
    else
        return j
    end
end

function mutate(rng::AbstractRNG, md::MutationDist, op::AbstractGeneOp)
    if rand(rng, md.d_mutate_op)
        rand(rng, md.op_inventory)
    else
        op
    end
end

function mutate(rng::AbstractRNG, md::MutationDist, inst::Instruction)
    op = mutate(rng, md, inst.op)
    rdd = RandomDuplicateDelete(
        rng,
        md.d_duplicate_index,
        md.d_delete_index,
        inst.operand_ixs)
    operand_ixs = mutate(rng, md, collect(rdd))
    return Instruction(op, operand_ixs)
end

function mutate(rng::AbstractRNG, md::MutationDist, v::AbstractVector)
    return map(x -> mutate(rng, md, x), v)
end

function mutate(rng::AbstractRNG, md::MutationDist, g::Genome)
    instruction_blocks = map(g.instruction_blocks) do block
        rdd = RandomDuplicateDelete(
            rng,
            md.d_duplicate_instruction,
            md.d_delete_instruction,
            block)
        mutate(rng, md, collect(rdd))
    end
    return Genome(instruction_blocks)
end
