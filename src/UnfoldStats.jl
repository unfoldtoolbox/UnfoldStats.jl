module UnfoldStats

using Unfold
using BSplineKit
using StatsModels
using LinearAlgebra
using Statistics # for mean

using StatsAPI
import StatsAPI: pvalue
include("pvalue.jl")
include("extract_coefs.jl")
include("cvMANOVA.jl")

# export functions to extract model coefficients 
export extract_coefs
export pvalue
export cvMANOVA

end
