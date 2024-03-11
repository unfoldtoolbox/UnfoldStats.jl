module UnfoldStats

using Unfold
using BSplineKit
using StatsModels

include("extract_coefs.jl")

# export functions to extract model coefficients 
export extract_coefs
end
