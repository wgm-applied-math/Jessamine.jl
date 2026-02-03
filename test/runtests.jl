include("TestGenome.jl")
include("TestRidgeLayer.jl")
include("TestEvolution.jl")
include("TestTimeSeries.jl")

function main()
    println("=== TestGenome ===")
    TestGenome.main()

    println("=== TestTimeSeries ===")
    TestTimeSeries.main()

    println("=== TestRidge ===")
    TestRidge.main()

    println("=== TestEvolution ===")
    TestEvolution.main()
end

main()
