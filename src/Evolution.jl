export SelectionSpec, SelectionDist, Population
export generation_size
export random_genome, random_initial_population, next_generation

"""
Parameters for tournament selection.
"""
@kwdef struct SelectionSpec
    "Keep this many high-rated agents from the current generation."
    num_to_keep::Int

    "Create this many new agents to add to the next generation."
    num_to_generate::Int

    """When choosing parents, pick two agents and random, and use
    the one with the highest rating with this probability.
    Otherwise, use the other one."""
    p_take_better::Float64
end

"""
    generation_size(s_spec::SelectionSpec)

Return the size of a generation.
"""
function generation_size(s_spec::SelectionSpec)
    return s_spec.num_to_keep + s_spec.num_to_generate
end

"""
Distribution object for tournament selection.
"""
struct SelectionDist
    spec::SelectionSpec
    d_take_highest::Bernoulli
end

function SelectionDist(s_spec::SelectionSpec)
    return SelectionDist(
        s_spec,
        Bernoulli(s_spec.p_take_better))
end

"""An agent has a rating, a genome, a parameter vector, and an
extra bit of data that depends on how rating is done."""
struct Agent{R,T}
    rating::R
    genome::Genome
    parameter::Vector
    extra::T
end

function short_show(io::IO, a::Agent)
    println(io, "rating = $(a.rating)")
    short_show(io, a.genome)
    println(io, "parameter = $(a.parameter)")
    println(io, "extra = $(a.extra)")
end

"""
    isless(x::Agent, y::Agent)

Compare the `rating` field of the two agents.
This definition allows a population to be sorted.
"""
Base.isless(x::Agent, y::Agent) = isless(x.rating, y.rating)

"""
A population is a single generation of living `Agent`s.
"""
struct Population
    agents::Vector{Agent}
end

function short_show(io::IO, p::Population)
    count = 0
    for j in eachindex(p.agents)
        a = p.agents[j]
        println(io, "slot $j:")
        short_show(io, a)
        println(io, "")
        count += 1
        if count >= 10
            println(io, "...")
            break
        end
    end
end

"""
    random_genome(rng::AbstractRNG, g_spec::GenomeSpec, m_spec::MutationSpec, arity_dist::Distribution)

Produce a random genome.  It will contain an instruction block
for each mutable slot in the state array as specified by
`g_spec`. Each block will contain one random instruction.  The
operator is picked uniformly at random from
`m_spec.op_inventory`.  The number of operands is drawn from
`arity_dist`.  The operands are drawn uniformly from the set of
possible indices.
"""
function random_genome(rng::AbstractRNG, g_spec::GenomeSpec,
        m_spec::MutationSpec, arity_dist::Distribution)
    index_max = workspace_size(g_spec)
    index_dist = DiscreteUniform(1, index_max)
    num_instruction_blocks = g_spec.output_size + g_spec.scratch_size
    instruction_blocks = Vector(undef, num_instruction_blocks)
    for j in 1:num_instruction_blocks
        op = rand(rng, m_spec.op_inventory)
        num_operands = rand(rng, arity_dist)
        operands = rand(rng, index_dist, num_operands)
        instruction_blocks[j] = [Instruction(op, operands)]
    end
    return Genome(instruction_blocks)
end

"""
    random_initial_population(rng::AbstractRNG, g_spec::GenomeSpec, m_spec::MutationSpec, arity_dist::Distribution, s_spec::SelectionSpec, grow_and_rate::Function)

Make a random initial population.  The number of genomes is
specified by adding the number of new genomes per generation and
the number of keepers given in `s_spec`.  Other parameters are
passed to `random_genome` to produce random genomes.
The `grow_and_rate` function is called with `rng, g_spec, genome`.
"""
function random_initial_population(
        rng::AbstractRNG,
        g_spec::GenomeSpec,
        m_spec::MutationSpec,
        arity_dist::Distribution,
        s_spec::SelectionSpec,
        grow_and_rate::Function)::Population
    pop_size = s_spec.num_to_keep + s_spec.num_to_generate
    agents = Vector(undef, pop_size)
    j = 1
    while j <= pop_size
        genome = random_genome(rng, g_spec, m_spec, arity_dist)
        agent = grow_and_rate(rng, g_spec, genome)
        if isnothing(agent)
            continue
        else
            agents[j] = agent
            j += 1
        end
    end
    return Population(agents)
end

"""
    pick_parent(rng::AbstractRNG, s_dist::SelectionDist, pop::Population)::Agent

Using tournament selection, pick a parent from the population.
"""
function pick_parent(rng::AbstractRNG, s_dist::SelectionDist, pop::Population)::Agent
    cs = rand(rng, pop.agents, 2)
    if rand(rng, s_dist.d_take_highest)
        return maximum(cs)
    else
        return minimum(cs)
    end
end

"""
    new_genome(rng::AbstractRNG, s_dist::SelectionDist, m_dist::MutationDist, pop::Population)

Produce a new offspring genome.
Pick two parents using tournament selection.
Recombine their genomes, and apply mutations.
"""
function new_genome(
        rng::AbstractRNG,
        s_dist::SelectionDist,
        m_dist::MutationDist,
        pop::Population)::Genome
    a1 = pick_parent(rng, s_dist, pop)
    a2 = pick_parent(rng, s_dist, pop)
    gr = recombine(rng, a1.genome, a2.genome)
    grm = mutate(rng, m_dist, gr)
    return grm
end

"""
    next_generation(rng::AbstractRNG, g_spec::GenomeSpec, s_dist::SelectionDist, m_dist::MutationDist, pop::Population, grow_and_rate::Function; sense=MinSense)

Produce the next generation of a population by selection and
mutation.  The offsprings' genomes are produced by `new_genome`,
and fed to `grow_and_rate(rng, genome)`, which should "grow" each
organism and return the rating and extra data as an `Agent`,
which is inserted into the population.  The `sense` parameter
specifies whether selection should minimize or maximize the
rating.
"""
function next_generation(
        rng::AbstractRNG,
        g_spec::GenomeSpec,
        s_dist::SelectionDist,
        m_dist::MutationDist,
        pop::Population,
        grow_and_rate::Function;
        sense = MinSense)::Population
    s_spec = s_dist.spec
    new_agents = Vector(undef, s_spec.num_to_generate)
    j = 1
    while j <= s_spec.num_to_generate
        agent = grow_and_rate(rng, g_spec, new_genome(rng, s_dist, m_dist, pop))
        if isnothing(agent)
            continue
        else
            new_agents[j] = agent
            j += 1
        end
    end
    keepers = pop.agents[1:(s_spec.num_to_keep)]
    next_rated_genomes = vcat(new_agents, keepers)
    rev = sense == Optimization.MaxSense
    sort!(next_rated_genomes, rev = rev)
    return Population(next_rated_genomes)
end
