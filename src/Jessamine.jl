module Jessamine

using Distributions
using LinearAlgebra
using Optim
using Optimization
using OptimizationOptimJL
using Random
using SciMLBase
using Symbolics
using Base.Threads

include("GenomeCore.jl")
include("RandomDuplicateDelete.jl")
include("Mutation.jl")
include("Recombination.jl")
include("Operations.jl")
include("Evolution.jl")
include("SymbolicForm.jl")
include("RidgeLayer.jl")

end
