using LinearAlgebra
using Optim
using Optimization
using OptimizationOptimJL
using Random

macro ignore(args...) end

include("TestGenome.jl")
include("TestRidgeLayer.jl")
include("TestEvolution.jl")
