using UnfoldStats
include("setup.jl")

@testset "UnfoldStats.jl" begin
    include("test-extract_coefs.jl")
end

@testset "MixedModelsPermutations+ClusterDepth" begin
    include("test-MixedModelsPermutations")
end