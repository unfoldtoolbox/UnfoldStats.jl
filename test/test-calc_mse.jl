@testset "MSE calculation" begin
    # Sample DataFrames
    results = DataFrame(eventname=["A", "A", "B", "B"], yhat=[1.0, 2.0, 3.0, 4.0], factor1=[1, 1, 2, 2], factor2=[1, 2, 1, 2])
    ground_truth = DataFrame(event=["A", "A", "B", "B"], yhat=[1.1, 1.9, 3.1, 3.9], factor1=[1, 1, 2, 2], factor2=[1, 2, 1, 2])
    effects_dict = Dict("factor1" => [1, 2], "factor2" => [1, 2])

    # Test with valid inputs
    @test calculate_mse(results, ground_truth, effects_dict) ≈ 0.01

    # Test with missing columns in results
    results_missing = DataFrame(eventname=["A", "A", "B", "B"], yhat=[1.0, 2.0, 3.0, 4.0], factor1=[1, 1, 2, 2])
    @test_throws AssertionError calculate_mse(results_missing, ground_truth, effects_dict)

    # Test with missing columns in ground_truth
    ground_truth_missing = DataFrame(event=["A", "A", "B", "B"], yhat=[1.1, 1.9, 3.1, 3.9], factor1=[1, 1, 2, 2])
    @test_throws AssertionError calculate_mse(results, ground_truth_missing, effects_dict)

    # Test with empty DataFrames
    empty_results = DataFrame()
    empty_ground_truth = DataFrame()
    @test_throws AssertionError calculate_mse(empty_results, empty_ground_truth, effects_dict)
end