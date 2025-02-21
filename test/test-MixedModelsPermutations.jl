


struct LMMClusterSimulation <: SimulationType end
PredictorStyle(::Type{LMMClusterSimulation}) = CategoricalPredictors()
ChannelStyle(::Type{LMMClusterSimulation}) = MultiChannel()
EventStyle(::Type{LMMClusterSimulation}) = SingleEventType()




β = [1, -2]
σs = Dict(:subject => [1, 0.1])#, 0.1])


sim_type = LMMClusterSimulation()
simulation = define_simulation(sim_type, β, σs, n_subjects = 5)
models = sim_and_fit(sim_type, simulation, UnfoldLinearMixedModel, seed = 1)
model = models.unfoldmodel[1]



ext = Base.get_extension(UnfoldStats, :UnfoldStatsMixedModelsPermutationsExt)
@testset "get_lmm_statistic" begin
    @test ext.get_lmm_statistic(model, 2, :β) ≈ last.(modelfit(model).β[2:2:end])
    @test length(ext.get_lmm_statistic(model, 2, :z)) ==
          length(ext.get_lmm_statistic(model, 2, :β))
end

@testset "lmm_permutations" begin

    permuted = ext.lmm_permutations(StableRNG(1), model, data, 2; n_permutations = 10)
    @test permuted ==
          ext.lmm_permutations(StableRNG(1), model, data, 2; n_permutations = 10) # test for repeated permutations
    @test size(permuted) == (3, 15, 10)

    _permuted = ext.lmm_permutations(
        StableRNG(1),
        model,
        data,
        2;
        n_permutations = 10,
        time_selection = 5:7,
    )
    @test size(_permuted) == (3, 3, 10)
    @test permuted[:, 5:7, :] == _permuted
end

@testset "lmm_clusterdepth_pvalues" begin
    permuted = ext.lmm_permutations(StableRNG(1), model, data, 2; n_permutations = 10)
    observed = ext.get_lmm_statistic(model, 2, :z)

    _observed = rand(length(observed)) .- 0.5
    @test all(
        ext.lmm_clusterdepth_pvalues(
            StableRNG(1),
            _observed,
            permuted;
            clusterforming_threshold = 1.8,
        ) .== 1.0,
    )

    _observed[4:6] .= 10
    @test all(
        ext.lmm_clusterdepth_pvalues(
            StableRNG(1),
            _observed,
            permuted;
            clusterforming_threshold = 1.8,
        )[4:6] .< 0.5,
    )

    # test clusterforming threshold
    @test all(
        ext.lmm_clusterdepth_pvalues(
            StableRNG(1),
            _observed,
            permuted;
            clusterforming_threshold = 10.8,
        ) .== 1.0,
    )
end

@testset "clusterdepth pvalue" begin

    p = pvalue(
        MersenneTwister(1),
        model,
        models.data[1],
        2;
        n_permutations = 20,
        clusterforming_threshold = 1.8,
    )
    @test size(p) == (3, 15)
    @test all(p[:, 4] .< 0.1)


    # pretty much nothing over the threshold for coefficient 1
    @test all(
        pvalue(
            MersenneTwister(1),
            model,
            models.data[1],
            1;
            n_permutations = 20,
            clusterforming_threshold = 1.8,
        ) .== 1.0,
    )

end