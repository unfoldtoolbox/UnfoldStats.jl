module UnfoldStats

using Unfold
using BSplineKit
using StatsModels, StatsBase
using MixedModelsSim # Only needed for factorproduct, maybe copy that function?

include("extract_coefs.jl")

# export functions to extract model coefficients 
export extract_coefs, calculate_mse
end
