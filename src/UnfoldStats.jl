module UnfoldStats

using BSplineKit: Reexport
using Unfold
using BSplineKit
using StatsModels
using StatsAPI
import StatsAPI: pvalue
include("pvalue.jl")
include("extract_coefs.jl")

# export functions to extract model coefficients 
export extract_coefs
export pvalue

end
