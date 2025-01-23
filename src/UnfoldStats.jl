module UnfoldStats

using BSplineKit: Reexport
using Unfold
using BSplineKit
using StatsModels
using StatsAPI
include("pvalue.jl")
include("extract_coefs.jl")

# export functions to extract model coefficients 
export extract_coefs

end
