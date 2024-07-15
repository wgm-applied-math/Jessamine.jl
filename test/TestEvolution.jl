println("Test next generation")

module EV
using Distributions
using LinearAlgebra
using Optimization
using Random
using Symbolics
using Test
using Jessamine
using ..RD

rng = Xoshiro(3208201)

g_spec = GenomeSpec(4, 0, 1, 2, 3)
index_max = workspace_size(g_spec)
m_spec = MutationSpec(
    p_mutate_op = 0.01,
    p_mutate_index = 0.1,
    p_duplicate_index = 0.0001,
    p_delete_index = 0.0001,
    p_duplicate_instruction = 0.00001,
    p_delete_instruction = 0.00001,
    op_inventory = [Add(), Multiply()])
m_dist = MutationDist(m_spec, index_max)
arity_dist = DiscreteNonParametric([1, 2, 3], [0.25, 0.5, 0.25])
s_spec = SelectionSpec(20, 30, 0.6)
gen_size = generation_size(s_spec)
s_dist = SelectionDist(s_spec)
lambda_b = 1e-6
lambda_p = 1e-6
lambda_o = 1e-6

# Adapt least_squares_ridge_grow_and_rate so that it can be used by next_generation
function grow_and_rate(rng, g_spec, genome)
    return least_squares_ridge_grow_and_rate(
        RD.xs, RD.y, lambda_b, lambda_p, lambda_o,
        g_spec, genome)
end

# Indices into the state vector
z1, z2, z3, z4, p1, x1, x2 = 1:index_max

# This should output [1, x1, x2, x1 * x2]
g_check = Genome(
    [[Instruction(Multiply(), Int[])], # creates a 1
    [Instruction(Add(), [x1])],
    [Instruction(Multiply(), [x2])],
    [Instruction(Multiply(), [x1, x2])]])

println("This should result in 2 + 3 x1 + 3 x2 - 3 x1 * x2")
a_check = grow_and_rate(rng, g_spec, g_check)
println("a_check =")
short_show(a_check)
println("end a_check")
println()
@show round.(coefficients(a_check.extra))
@test round.(coefficients(a_check.extra)) == [2, 3, 3, -3]
@test intercept(a_check.extra) == 0

pop_init = random_initial_population(rng, g_spec, m_dist, arity_dist, s_spec, grow_and_rate)
println("pop_init = ")
short_show(pop_init)
println("end pop_init")
println("")

function ev_loop(pop_init::Population, num_generations::Int)
    pop_next = pop_init
    for t in 1:num_generations
        pop_next = next_generation(
            rng, g_spec, s_dist, m_dist, pop_next, grow_and_rate; sense = MinSense)
        if mod(t, 10) == 0
            println("Generation $t, best = $(pop_next.agents[1].rating)")
        end
    end
    return pop_next
end

@time @eval pop_next = ev_loop(pop_init, 100)
println("pop_next = ")
short_show(pop_next)
println("end pop_next")
println("")

println("Symbolic form of best agent")

a_best = pop_next.agents[1]

sym_res = model_symbolic_output(g_spec, a_best)
x = sym_res.x

short_show(a_best)
@show sym_res.y_sym
@show a_best.parameter
@show a_best.extra
@show sym_res.y_num

y_num = simplify(sym_res.y_num; expand = true)
@show y_num

c1 = Symbolics.coeff(y_num)
cx1 = Symbolics.coeff(Symbolics.coeff(y_num, x[1]))
cx2 = Symbolics.coeff(Symbolics.coeff(y_num, x[2]))
cx1x2 = Symbolics.coeff(Symbolics.coeff(y_num, x[2]), x[1])
cvec = [c1, cx1, cx2, cx1x2]
@show cvec
@show round.(cvec)
@test round.(cvec) == [2, 3, 3, -3]
end
