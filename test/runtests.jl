using UnfoldStats
include("setup.jl")

@testset "UnfoldStats.jl" begin
    include("extract_coefs.jl")
end

include("test-calc_mse.jl")