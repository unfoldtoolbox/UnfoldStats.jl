using StableRNGs
using LinearAlgebra
using Unfold
using UnfoldSim
@testset "cvMANOVA" begin
    rng = StableRNG(1)

    @testset "helpers" begin
        Xi = randn(rng, 20, 3)
        Sinv = UnfoldStats.calculate_Σinv(Xi; λ = 0.1)
        @test size(Sinv) == (3, 3)
        @test Sinv ≈ Sinv'

        CCXXCC = Matrix{Float64}(I, 2, 2)
        Btrain = reshape([1.0, 2.0, 3.0, 4.0], 2, 2)
        Btest = reshape([2.0, 1.0, 0.0, 1.0], 2, 2)
        d = UnfoldStats.cvmanova_D(CCXXCC, Matrix{Float64}(I, 2, 2), Btrain, Btest)
        @test d isa Real
        @test d ≈ 8.0
    end

    @testset "array API" begin
        n_train = 24
        n_test = 12
        n_chan = 3
        n_time = 8
        n_coef = 2

        X_train = randn(rng, n_train, n_coef)
        X_test = randn(rng, n_test, n_coef)
        β_train = randn(rng, n_chan, n_time, n_coef)
        β_test = randn(rng, n_chan, n_time, n_coef)
        Y_train = randn(rng, n_chan, 5, n_train)

        C = [-1, 1] #Matrix{Float64}(I, n_coef, n_coef)
        D = UnfoldStats.cvMANOVA(
            β_train,
            β_test,
            X_train,
            X_test,
            Y_train;
            C = C,
            Σinv_kwargs = (; λ = 0.05),
        )

        @test D isa AbstractArray
        @test length(D) == n_time
        @test all(isfinite, D)
    end

    @testset "model API" begin
        eeg, evts = UnfoldSim.predef_eeg(; return_epoched = true, multichannel = true)
        f = @formula 0 ~ 0 + condition

        solver = Unfold.solver_cv(n_folds = 3, shuffle = false, fit_test = true)

        m = fit(
            UnfoldModel,
            f,
            evts,
            eeg,
            1:size(eeg, 2);
            solver = solver,
            contrasts = Dict(:condition => Unfold.StatsModels.FullDummyCoding()),
        )


        D_all = cvMANOVA(m, eeg; C = [-1, 1])
        @test length(D_all) == length(modelfit(m).folds[1])

        D_first = cvMANOVA(m, eeg, 1; C = [-1, 1])
        @test length(D_first) == size(modelfit(m).estimate, 2)

    end
end
