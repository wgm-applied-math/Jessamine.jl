module Jessamine

using Dates
using Distributions
using LinearAlgebra
using MLJ
using MLJLinearModels
using MLJModelInterface
using NamedTupleTools
using Optim
using Optimization
using OptimizationOptimJL
using PyCall
using Random
using SciMLBase
using Symbolics
using SymPy
using SymPyCore
using Tables
using TermInterface
using Base.Threads

include("TablesUtil.jl")
include("GenomeCore.jl")
include("RandomDuplicateDelete.jl")
include("Mutation.jl")
include("Recombination.jl")
include("Operations.jl")
include("Evolution.jl")
include("SymbolicForm.jl")
include("SymPyForm.jl")
include("AbstractModelLayer.jl")
include("AbstractLinearModelLayer.jl")
include("RidgeLayer.jl")
include("MachineLayer.jl")
include("MLJInside.jl")
include("MLJOutside.jl")
include("TimeSeriesLayer.jl")
include("TimeSeriesRidgeLayer.jl")

end
