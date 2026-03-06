module TestRidge
using LinearAlgebra
using Test
using Jessamine

include("RandomData.jl")

z1, z2, z3, t1, t2, p1, pm1, x1, x2 = 1:9

genome = Genome(
    [[Instruction(Add(), [x2])],
     [Instruction(Multiply(), [x1, t2])],
     [Instruction(Multiply(), [])],
     [Instruction(Multiply(), [x2, pm1])],
     [Instruction(Add(), [p1, t1])],
]
)

parameter = [1.0, -1.0]

g_spec = GenomeSpec(3, 2, 2, 2, 4)

g_res = run_genome_to_last(g_spec, genome, parameter, RD.xs)

function main()
    @show least_squares_ridge(RD.xs, RD.y, 0.0, g_spec, genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 1e-6, g_spec, genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 0.1, g_spec, genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 0.5, g_spec, genome, parameter)
    @show least_squares_ridge(RD.xs, RD.y, 1.0, g_spec, genome, parameter)

    global lambda_b = 1e-6
    global lambda_p = lambda_b
    global lambda_o = lambda_b
    global agent = least_squares_ridge_grow_and_rate(
        RD.xs, RD.y, lambda_b, lambda_p, lambda_o, g_spec, genome)
    global y_rgr_pred = model_predict(g_spec, agent, RD.xs)
    global residuals = RD.y - y_rgr_pred

    @show agent.parameter
    @show agent.extra
    @show y_rgr_pred
end

end # module TestRidge
