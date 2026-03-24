using UnfoldStats
include("setup.jl")

@testset "UnfoldStats.jl" begin
    include("test-extract_coefs.jl")

end
@testset "cvMANOVA" begin
    include("test-cvMANOVA.jl")
end
@testset "MixedModelsPermutations+ClusterDepth" begin
    include("test-MixedModelsPermutations.jl")
end
