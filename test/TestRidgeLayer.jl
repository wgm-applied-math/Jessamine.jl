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
using LinearAlgebra
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

spec = GenomeSpec(2, 3, 2, 2, 4)

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
    y_rgr_pred = linear_model_predict(RG.spec, agent, RD.xs)
    residuals = RD.y - y_rgr_pred
    return (agent, y_rgr_pred, residuals)
end

agent, y_rgr_pred, rgr_residuals = test_ridge_grow_and_rate()

@show agent.parameter
@show agent.extra
@show y_rgr_pred

end # module RG

println("Test symbolics")

module RS
using LinearAlgebra
using Symbolics
using Jessamine
using ..RD
using ..RG

# Things work best if we use an array of Symbolic objects
# instead of a Symbolic array object.
x = Symbolics.variables(:x, 1:2)

z_sym = run_genome(RG.spec, RG.genome, Num.(RG.agent.parameter), x)[end]
b = coefficients(RG.agent.extra)
y_pred_sym = dot(z_sym, b)
@show y_pred_sym
y_pred_simp = simplify(y_pred_sym)
@show y_pred_simp

end # module RS

println("End Test ridge regression layer")
println("----------")
