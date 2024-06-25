println("Test ridge regression layer")

module RD
using Distributions
using Random

rng = Xoshiro(161004)

num_points = 30
x1_dist = Normal(0.0, 1.0)
x1 = rand(rng, x1_dist, num_points)
x2_dist = Normal(0.0, 2.0)
x2 = rand(rng, x2_dist, num_points)
y = @. 2 + 3 * (x2 + x1 * (1 - x2))
# which is 2 + 3 x2 + 3 x1 ( 1 - x2)
#          2 + 3 x2 + 3 x1 - 3 x1 x2

xs = [x1, x2]
end # module RD

module RG
using Test
using Jessamine
using ..RD

y1, t1, t2, t3, p1, pm1, x1, x2 = 1:8

genome = Genome(
    [[Instruction(Add(), [x2, t3])],
    [Instruction(Multiply(), [pm1, x2])],
    [Instruction(Add(), [p1, t1])],
    [Instruction(Multiply(), [x1, t2])]
]
)

parameter = [1.0, -1.0]

spec = GenomeSpec(1, 3, 2, 2, 4)

g_res = run_genome(spec, genome, parameter, RD.xs)
y = g_res[end][1]

function test_ridge()
    @show least_squares_ridge(RD.xs, RD.y, 0.0, RG.spec, RG.genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 1e-6, RG.spec, RG.genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 0.1, RG.spec, RG.genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 0.5, RG.spec, RG.genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 1.0, RG.spec, RG.genome, parameter)
end

test_ridge()

function test_ridge_grow_and_rate()
    lambda_b = 1e-6
    lambda_p = lambda_b
    lambda_o = lambda_b
    agent = least_squares_ridge_grow_and_rate(
        RD.xs, RD.y, lambda_b, lambda_p, lambda_o, RG.spec, RG.genome)
    g_res = run_genome(RG.spec, RG.genome, agent.parameter, RD.xs)
    y_rgr = g_res[end][1]
    b = agent.extra
    y_rgr_pred = @. b[1] + b[2] * y_rgr
    residuals = RD.y - y_rgr_pred
    return (agent, y_rgr_pred, residuals)
end

agent, y_rgr_pred, rgr_residuals = test_ridge_grow_and_rate()

@show agent.parameter
@show agent.extra
@show y_rgr_pred

@test round(agent.extra[1]) == 2
@test round(agent.extra[2]) == 3

end # module RG

println("Test symbolics")

module RS
using LinearAlgebra
using Symbolics
using Jessamine
using ..RD
using ..RG

# The documentation says that fewer things work with symbolic
# arrays like this, as opposed to an array of symbols.  Not sure
# if that's importat yet.
@variables x[1:2]

y_rgr_sym = run_genome(RG.spec, RG.genome, RG.agent.parameter, x)[end][1]
b = RG.agent.extra
y_pred_sym = b[1] + dot(b[2:end], y_rgr_sym)
@show y_pred_sym
y_pred_simp = simplify(y_pred_sym)
@show y_pred_simp

end # module RS

println("End Test ridge regression layer")
println("----------")
