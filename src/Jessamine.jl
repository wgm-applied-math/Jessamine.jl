module Jessamine

using DataFrames
using Distributions
using LinearAlgebra
using MLJ
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
#include("MachineLayer.jl")
#include("MLJInside.jl")

end
