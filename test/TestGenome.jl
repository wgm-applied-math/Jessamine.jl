module T
using Distributions
using Random
using Test

using Jessamine

rng = Xoshiro(202406181224)

println("Test run_genome with the Henon map")

a = 1.4
ma = -a
b = 0.3
input = [0.5, 0.7]

f((x, y)) = @. [1 - a * x^2 + y, b * x]

function g((x, y))
    t1 = @. ma * x * x
    x1 = @. 1 + t1 + y
    y1 = @. b * x
    return [x1, y1]
end

@test f(input) == g(input)

g_spec = GenomeSpec(2, 1, 3, 2, 2)

# Map from variables to indices
x1, y1, t1, p1, pma, pb, x0, y0 = 1:8

instruction_blocks = [
    [Instruction(Add(), [p1, t1, y0])],
    [Instruction(Multiply(), [pb, x0])],
    [Instruction(Multiply(), [pma, x0, x0])]
]

parameter = [1.0, ma, b]

genome = Genome(instruction_blocks)

function h(input)
    output = run_genome(g_spec, genome, parameter, input)
    return output
end

@show h(input)

@test f(input) == h(input)[end]

println("Test syntax")
@show genome_as_syn = to_syntax(g_spec, genome)
g_comp = Jessamine.eval(genome_as_syn)
function hc(input)
    output = run_compiled_genome(g_spec, g_comp, parameter, input)
    return output
end
@test f(input) == hc(input)[end]

println("Test vectorization")

function hv(input)
    output = run_genome(g_spec, genome, parameter, input)
    return output
end

function hvc(input)
    output = run_compiled_genome(g_spec, g_comp, parameter, input)
    return output
end

v_input = [[0.5, 0.5, 0.2], [0.7, 0.8, 0.8]]

println("Test vectorization:")
@show hv(v_input)
@test f(v_input) == hv(v_input)[end]
@test f(v_input) == hvc(v_input)[end]

println("Test RandomDuplicateDelete")

rdd_empty = RandomDuplicateDelete(rng, Bernoulli(0.2), Bernoulli(0.2), [])
@test nothing === iterate(rdd_empty)

# Test empty iteration, take 1
for t in rdd_empty
    @test false
end

# Test empty iteration, take 2
@test [] == [x for x in rdd_empty]

rdd = RandomDuplicateDelete(rng, Bernoulli(0.2), Bernoulli(0.2), 1:100)
rdd_res = [x for x in rdd]
@show rdd_res
@show collect(rdd)

println("----------")

println("Test mutation")

m_spec = MutationSpec(
    p_mutate_op = 0.5,
    p_mutate_index = 0.5,
    p_duplicate_index = 0.2,
    p_delete_index = 0.2,
    p_duplicate_instruction = 0.2,
    p_delete_instruction = 0.2,
    op_inventory = [Add(), Multiply()])

m_dist = MutationDist(m_spec, 9)

inst1 = Instruction(Add(), [1, 3, 5])
@show inst1
inst2 = mutate(rng, m_dist, inst1)
@show inst2
inst3 = mutate(rng, m_dist, inst2)
@show inst3
inst4 = mutate(rng, m_dist, inst3)
@show inst4
inst5 = mutate(rng, m_dist, inst4)
@show inst5
inst6 = mutate(rng, m_dist, inst5)
@show inst6
inst7 = mutate(rng, m_dist, inst6)
@show inst7
inst8 = mutate(rng, m_dist, inst7)
@show inst8
inst9 = mutate(rng, m_dist, inst8)
@show inst9
inst10 = mutate(rng, m_dist, inst9)
@show inst10

genome1 = Genome(
    [
    [inst1, inst2],
    [inst3, inst4, inst5],
    [inst6],
    [],
    [inst7, inst8, inst9, inst10]])
println("genome1 = ")
short_show(genome1)
genome2 = mutate(rng, m_dist, genome1)
println("genome2 = ")
short_show(genome2)

println("----------")

println("Test recombination")

genome3 = Genome(
    [
    [inst9, inst10, inst1],
    [inst6, inst7, inst8],
    [inst6],
    [inst1, inst2],
    [inst2, inst3, inst4]])
println("genome3 = ")
short_show(genome3)
genome4 = recombine(rng, genome1, genome3)
println("genome4 = ")
short_show(genome4)
genome5 = recombine(rng, genome3, genome1)
println("genome5 = ")
short_show(genome5)

println("----------")

end # module T
