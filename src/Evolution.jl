# Make language server happy
if false
    using Dates
    using Distributions
    include("GenomeCore.jl")
    include("Operations.jl")
    include("Recombination.jl")
    include("Mutation.jl")
end

export SelectionSpec, SelectionDist, Population
export Agent
export generation_size
export random_genome, random_initial_population, next_generation
export EvolutionSpec, evolution_loop
export vns_evolution_loop

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
    p_take_better::Float64 = 0.6

    """When choosing parents, pick the top-rated individual
    with this probability first, otherwise proceed to tournament
    selection."""
    p_take_very_best::Float64 = 0.0
end

function Base.convert(::Type{SelectionSpec}, s_spec::SelectionSpec)
    return s_spec
end

function Base.convert(::Type{SelectionSpec}, x)
    return SelectionSpec(
    x.num_to_keep,
    x.num_to_generate,
    x.p_take_better,
    x.p_take_very_best)
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
    d_take_better::Bernoulli
    d_take_very_best::Bernoulli
end

function SelectionDist(s_spec)
    return SelectionDist(
        convert(SelectionSpec, s_spec),
        Bernoulli(s_spec.p_take_better),
        Bernoulli(s_spec.p_take_very_best))
end

"""An agent has a rating, a genome, a parameter vector, and an
extra bit of data that depends on how rating is done."""
struct Agent{R, G <: AbstractGenome, VPar <: AbstractVector, T}
    rating::R
    genome::G
    parameter::VPar
    extra::T
end

function compile(g_spec, agent)
    return Agent(agent.rating, compile(g_spec, agent.genome), agent.parameter, agent.extra)
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
        m_dist::MutationDist, arity_dist::Distribution)
    index_max = workspace_size(g_spec)
    index_dist = DiscreteUniform(1, index_max)
    num_instruction_blocks = g_spec.output_size + g_spec.scratch_size
    instruction_blocks = Vector(undef, num_instruction_blocks)
    for j in 1:num_instruction_blocks
        op = m_dist.op_inventory[rand(rng, m_dist.d_op)]
        num_operands = rand(rng, arity_dist)
        operands = rand(rng, index_dist, num_operands)
        instruction_blocks[j] = [Instruction(op, operands)]
    end
    g = Genome(instruction_blocks)
    # cg = compile(g_spec, g)
    return g
end

"""
Parameters needed to run `random_initial_population` and `evolution_loop`.
"""
@kwdef struct EvolutionSpec
    g_spec::GenomeSpec
    m_dist::MutationDist
    s_dist::SelectionDist
    grow_and_rate::Any
    num_generations::Int
    stop_on_innovation::Bool = false
end

function EvolutionSpec(
        g_spec::GenomeSpec,
        m_spec::MutationSpec,
        s_spec::SelectionSpec,
        grow_and_rate,
        num_generations::Int,
        stop_on_innovation::Bool = false
)
    src_index_max = workspace_size(g_spec)
    m_dist = MutationDist(m_spec, src_index_max)
    s_dist = SelectionDist(s_spec)
    return EvolutionSpec(
        g_spec,
        m_dist,
        s_dist,
        grow_and_rate,
        num_generations,
        stop_on_innovation
    )
end

function is_better(reference, other, sense)
    if sense == MinSense
        return other < reference
    else
        return other > reference
    end
end

"""
    random_initial_population(rng::AbstractRNG, g_spec::GenomeSpec, m_dist::MutationDist, arity_dist::Distribution, s_spec::SelectionSpec, grow_and_rate)

Make a random initial population.  The number of genomes is
specified by adding the number of new genomes per generation and
the number of keepers given in `s_spec`.  Other parameters are
passed to `random_genome` to produce random genomes.
The `grow_and_rate` function is called with `rng, g_spec, genome`.
"""
function random_initial_population(
        rng::AbstractRNG,
        g_spec::GenomeSpec,
        m_dist::MutationDist,
        arity_dist::Distribution,
        s_spec::SelectionSpec,
        grow_and_rate;
        sense = MinSense)::Population
    pop_size = s_spec.num_to_keep + s_spec.num_to_generate
    agents = Vector(undef, pop_size)
    Threads.@threads for j in eachindex(agents)
        agent = nothing
        while isnothing(agent)
            genome = random_genome(rng, g_spec, m_dist, arity_dist)
            # @assert check_genome(genome) "random init"
            agent = grow_and_rate(rng, g_spec, genome)
        end
        agents[j] = agent
    end
    rev = sense == Optimization.MaxSense
    sort!(agents, rev = rev)
    return Population(agents)
end

function random_initial_population(
        rng::AbstractRNG,
        e_spec::EvolutionSpec,
        arity_dist::Distribution,
        grow_and_rate;
        sense = MinSense)::Population

    return random_initial_population(
    rng,
    e_spec.g_spec,
    e_spec.m_dist,
    arity_dist,
    e_spec.s_dist.spec,
    grow_and_rate;
    sense = sense)
end

"""
    pick_parent(rng::AbstractRNG, s_dist::SelectionDist, pop::Population)::Agent

Using tournament selection, pick a parent from the population.
"""
function pick_parent(rng::AbstractRNG, s_dist::SelectionDist, pop::Population)::Agent
    if s_dist.spec.p_take_very_best > 0 && rand(rng, s_dist.d_take_very_best)
        return pop.agents[1]
    else
        cs = rand(rng, pop.agents, 2)
        if rand(rng, s_dist.d_take_better)
            return maximum(cs)
        else
            return minimum(cs)
        end
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
        pop::Population
)::Genome
    a1 = pick_parent(rng, s_dist, pop)
    a2 = pick_parent(rng, s_dist, pop)
    gr = recombine(rng, a1.genome, a2.genome)
    # @assert check_genome(gr) "gr"
    grm = mutate(rng, m_dist, gr)
    # @assert check_genome(grm) "grm"
    return grm
end

"""
    next_generation(rng::AbstractRNG, g_spec::GenomeSpec, s_dist::SelectionDist, m_dist::MutationDist, pop::Population, grow_and_rate::Function; sense=MinSense)

Produce the next generation of a population by selection and
mutation.  The offsprings' genomes are produced by `new_genome`,
and fed to `grow_and_rate(rng, g_spec, genome)`, which should "grow" each
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
        grow_and_rate;
        sense = MinSense)::Population
    s_spec = s_dist.spec
    new_agents = Vector(undef, s_spec.num_to_generate)
    Threads.@threads for j in eachindex(new_agents)
        agent = nothing
        while isnothing(agent)
            ng = new_genome(rng, s_dist, m_dist, pop)
            # cg = compile(g_spec, g)
            agent = grow_and_rate(rng, g_spec, ng)
        end
        new_agents[j] = agent
    end
    keepers = pop.agents[1:(s_spec.num_to_keep)]
    next_rated_genomes = vcat(new_agents, keepers)
    rev = sense == Optimization.MaxSense
    sort!(next_rated_genomes, rev = rev)
    return Population(next_rated_genomes)
end


"""
    evolution_loop(rng::AbstractRNG, g_spec::GenomeSpec, m_dist::MutationDist, s_dist::SelectionDist, grow_and_rate, num_generations::Integer, sense, pop_init::Population; stop_threshold=nothing, stop_on_innovation=false, generation_mod=10, enable_logging=true)

Advance `pop_init` by `num_generations` and return the new `Population`.

The various spec records and `grow_and_rate` function are passed
to `next_generation()`.

`sense` should be `MinSense` or `MaxSense` to
indicate whether the search should minimize or maximize the
rating.

If `stop_on_innovation` is `true`, the process stops as soon as
an agent with a better rating than the previous best is
discovered.

Progress is reported as an `@info` log message when the
generation number is a multiple of `generation_mod`.  If this
number is not > 0, this progress report is disabled.

Progress is reported as an `@info` log message when an agent with
a better rating than the previous best is discovered.

Setting `enable_logging` to `false` disables both the periodic
progress report and innovation messages.

If `stop_threshold` is a number rather than `nothing`, the
process stops as soon as an agent with a rating better than
`stop_threshold` is discovered.
"""

function evolution_loop(
        rng::AbstractRNG,
        g_spec::GenomeSpec,
        m_dist::MutationDist,
        s_dist::SelectionDist,
        grow_and_rate,
        num_generations::Integer,
        pop_init::Population;
        sense = MinSense,
        stop_threshold = nothing,
        stop_on_innovation = false,
        generation_mod = 10,
        enable_logging = true
)
    pop_next = pop_init
    best_rating = pop_init.agents[1].rating
    for t in 1:num_generations
        if !isnothing(stop_threshold) && is_better(stop_threshold, best_rating, sense)
            break
        end
        pop_next = next_generation(
            rng, g_spec, s_dist, m_dist, pop_next, grow_and_rate)
        best_in_gen = pop_next.agents[1]
        if is_better(best_rating, best_in_gen.rating, sense)
            best_rating = best_in_gen.rating
            if enable_logging
                @info "$(now()): Generation $t, new best = $best_rating"
            end
            if stop_on_innovation
                break
            end
        end
        if enable_logging && t > 0 && mod(t, generation_mod) == 0
            @info "$(now()): Generation $t, best = $best_rating"
        end
    end
    return pop_next
end

"""
    evolution_loop(
    rng::AbstractRNG,
    e_spec::EvolutionSpec
    pop_init::Population;
    sense=MinSense,
    stop_threshold=nothing,
    generation_mod=10,
    enable_logging=true
    )

Advance `pop_init` and return the new `Population`.

Calls `evolution_loop()` with the specifications unpacked from `e_spec`.

Progress is reported as an `@info` log message when the
generation number is a multiple of `generation_mod`.  If this
number is not > 0, this progress report is disabled.

Progress is reported as an `@info` log message when an agent with
a better rating than the previous best is discovered.

Setting `enable_logging` to `false` disables both the periodic
progress report and innovation messages.
"""

function evolution_loop(
        rng::AbstractRNG,
        e_spec::EvolutionSpec,
        pop_init::Population;
        sense = MinSense,
        stop_threshold = nothing,
        generation_mod = 10,
        enable_logging = true
)
    return evolution_loop(
        rng,
        e_spec.g_spec,
        e_spec.m_dist,
        e_spec.s_dist,
        e_spec.grow_and_rate,
        e_spec.num_generations,
        pop_init,
        sense = sense,
        stop_threshold = stop_threshold,
        stop_on_innovation = e_spec.stop_on_innovation,
        generation_mod = generation_mod,
        enable_logging = enable_logging
    )
end

"""
    vns_evolution_loop(
    rng::AbstractRNG,
    neighborhoods::AbstractVector{EvolutionSpec},
    num_epochs::Integer,
    pop_init::Population;
    sense=MinSense,
    stop_threshold=nothing,
    generation_mod=10,
    enable_logging=true
    )

Perform up to `num_epochs` iterations of variable neighborhood search.

Evolution is performed using `neighborhoods[1]` using `evolution_loop`.
If a new best rating is achieved, continue.
Otherwise, switch to using `neighborhoods[2]` etc.
Once a new best rating is achieved, return to `neighborhoods[1]`.

The keyword arguments are passed to `evolution_loop`.

If `enable_logging` is `true`, `@info` log messages are
produced after each epoch completes.
"""

function vns_evolution_loop(
        rng::AbstractRNG,
        neighborhoods::AbstractVector{EvolutionSpec},
        num_epochs::Integer,
        pop_init::Population;
        sense = MinSense,
        stop_threshold = nothing,
        generation_mod = 10,
        enable_logging = true
)
    pop_next = pop_init
    best_rating = pop_init.agents[1].rating
    neighborhood_index = 1
    for e in 1:num_epochs
        if !isnothing(stop_threshold) && is_better(stop_threshold, best_rating, sense)
            break
        end
        pop_next = evolution_loop(
            rng,
            neighborhoods[neighborhood_index],
            pop_next;
            sense = sense,
            stop_threshold = stop_threshold,
            generation_mod = generation_mod,
            enable_logging = enable_logging
        )
        best_in_epoch = pop_next.agents[1]
        if is_better(best_rating, best_in_epoch.rating, sense)
            best_rating = best_in_epoch.rating
            # Return to neighborhood 1
            neighborhood_index = 1
            if enable_logging
                @info "$(now()): Epoch $e, new best = $best_rating, returning to first neighborhood"
            end
        elseif neighborhood_index < length(neighborhoods)
            # Advance to the next neighborhood
            neighborhood_index += 1
            if enable_logging
                @info "$(now()): Epoch $e completed with no new best rating, advancing to neighborhood $neighborhood_index"
            end
        else
            if enable_logging
                @info "$(now()): Epoch $e completed with no new best rating, continuing in neighborhood $neighborhood_index"
            end
        end
    end
    return pop_next
end
