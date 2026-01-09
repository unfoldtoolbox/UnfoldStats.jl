@testset "extract_coefs" begin
    #---
    sim_type = UnitTestSimulation()

    # Print termnames to check how many β and σs values need to be defined
    # print(StatsModels.termnames(get_formula(sim_type).rhs))

    β = [1, 0, 0, 0, 0, 0]
    σs = Dict(:subject => [1, 0.1, 0.1, 0.1])

    simulation = define_simulation(sim_type, β, σs, n_subjects = 5)
    for model_type in [UnfoldLinearModel, UnfoldLinearModelContinuousTime]

        models = sim_and_fit(sim_type, simulation, model_type, seed = 1)
        model_1 = models[1, :unfoldmodel]

        symbols = UnfoldStats.extract_symbol(Unfold.formulas(model_1))

        if model_type == UnfoldLinearModel

            # Test helper functions
            @test symbols == [
                "(Intercept)",
                :continuous,
                "(Intercept)",
                :continuous,
                :continuous,
                :continuous,
                :continuous,
                :condition,
                :pet,
                (:condition, :pet),
            ]
            @test UnfoldStats.contained_or_equal.(:continuous, symbols) ==
                  [0, 1, 0, 1, 1, 1, 1, 0, 0, 0]

            @test UnfoldStats.get_predictor_string(:continuous) == ":continuous"
            @test UnfoldStats.get_predictor_string("(Intercept)") == "\"(Intercept)\""
            @test UnfoldStats.get_predictor_string((:condition, :pet)) ==
                  "(:condition, :pet)"

            @test Unfold.get_basis_names(model_1) ==
                  reduce(vcat, fill.(["stim", "fix"], [2, 8]))

            # Test exceptions
            @test_throws ArgumentError extract_coefs(model_1, "continuous", "fix")
            @test_throws ArgumentError extract_coefs(
                model_1,
                :continuous,
                "this-basis-doesnt-exist",
            )
            @test_throws ArgumentError extract_coefs(model_1, :pet, "stim")

            # Test `extract_coefs` method for single Unfold model
            coefs_1 = extract_coefs(model_1, :continuous, "fix")
            @test size(coefs_1) == (3, 14, 4)
            @test coefs_1[:] ==
                  subset(
                coeftable(model_1),
                :eventname => ByRow(==("fix")),
                :coefname => ByRow(
                    x -> x in [
                        "continuous",
                        "spl(continuous,1)",
                        "spl(continuous,2)",
                        "spl(continuous,3)",
                    ],
                ),
            ).estimate

            # Test `extract_coefs` method for a vector of Unfold models
            coefs = extract_coefs(models.unfoldmodel, :continuous, "stim")
            @test coefs[:, :, 1, 4][:] ==
                  subset(
                coeftable(models[4, :unfoldmodel]),
                :eventname => ByRow(==("stim")),
                :coefname => ByRow(==("continuous")),
            ).estimate

        else
            # Test helper functions
            @test length(symbols) == 1210
            @test Unfold.get_basis_names(model_1) |>
                  x ->
                vcat(x...) == reduce(vcat, fill.(["stim", "fix"], [161 * 2, 111 * 8]))

            # Test `extract_coefs` method for single Unfold model
            coefs_1 = extract_coefs(model_1, :continuous, "fix")
            @test size(coefs_1) == (3, 111, 4)
            @test coefs_1[:, 1, :][:] ==
                  subset(
                coeftable(model_1),
                :eventname => ByRow(==("fix")),
                :coefname => ByRow(
                    x -> x in [
                        "continuous",
                        "spl(continuous,1)",
                        "spl(continuous,2)",
                        "spl(continuous,3)",
                    ],
                ),
                :time => (ByRow(==(-0.1))),
            ).estimate

            # Test `extract_coefs` method for a vector of Unfold models
            coefs = extract_coefs(models.unfoldmodel, :continuous, "stim")
            @test coefs[:, :, 1, 4][:] ==
                  subset(
                coeftable(models[4, :unfoldmodel]),
                :eventname => ByRow(==("stim")),
                :coefname => ByRow(==("continuous")),
            ).estimate
        end
    end
end
