"""
    calculate_mse(results:DataFrame, ground_truth::DataFrame, effects_dict::Dict; weighted_mse::Bool)

Calculates the (weighted) mean squared error (MSE) of marginalized effects and their ground truth.

The function is specific for calculating MSE values between results obtained from Unfold.jl and their ground truth (EffectsDesign) simulated using UnfoldSim.
Note that currently, averages are weighted by covariates per event. In the future this is likely to change to weight values by size(covariate), e.g. spline predictors are weighted by the number of splines.

# Arguments (if needed)
- `results::DataFrame`: Results DataFrame obtained from `Unfold.effects()`.
- `ground_truth::DataFrame`: EffectsDesign DataFrame obtained from `UnfoldSim.effects()`.

# Keywordarguments
- `weighted_mse::Bool`: If true (default) then the function returns a weighted average, otherwise a vector containing one MSE value per event.

# Returns
if weighted_mse == true
- `overall_mse` : Mean squared error.
if weighted_mse == false
- `event_mse` : Vector of average MSE per event_mse.
- `event_weights` : weights to calculated weighted mse.
"""

function calculate_mse(results::DataFrame, ground_truth::DataFrame, effects_dict::Dict; weighted_mse:Bool=true)
    # Ensure the dataframes have effects_dict entries
    if typeof(first(keys(effects_dict))) != String
        tmp = tostringdict(effects_dict)
        @assert all(x -> x in names(results), collect(keys(tmp))) "Effects not found in results "
        @assert all(x -> x in names(ground_truth), collect(keys(tmp))) "Effects not found in ground_truth "
    else
        @assert all(x -> x in names(results), collect(keys(effects_dict))) "Effects not found in results "
        @assert all(x -> x in names(ground_truth), collect(keys(effects_dict))) "Effects not found in ground_truth "
    end

    # Add event collumn to make dfs consistent
    results.event = results.eventname

    # Add event to dict to make dict_list
    full_factorial = factorproduct(((; k => v) for (k, v) in pairs(effects_dict))...) |> DataFrame
    @debug "Full Factorial DataFrame:" full_factorial
    
    # In singular events ground_truth doesn't have an event collumn
    if !("event" ∈ names(ground_truth))
        ground_truth[:, :event] .= Any
    end

    # Group DataFrames by event
    group_results = groupby(results, :event)
    group_ground_truth = groupby(ground_truth, :event)
    @debug first(group_results, 2)
    @debug first(full_factorial, 2)
    event_mse = []
    event_weights::AbstractVector{<:Real} = Float64[]
    for both_groups in zip(group_results, group_ground_truth)
        tmp_mse = []
        for row in eachrow(full_factorial)
            collums = Symbol.(names(row))
            tmp_value = values(row)
            
            @debug "Current Factorial Row:" row
            @debug "Current Factorial Columns:" collums
            @debug "Current Factorial Values:" tmp_value

            #error("This is a debug message to check the current factorial row, columns, and values.")
            # Filter both groups and calculate MSE on current factorial 
            d1 = subset(both_groups[1], [col => x -> x .== val for (col, val) in zip(collums, tmp_value)])
            d2 = subset(both_groups[2], [col => x -> x .== val for (col, val) in zip(collums, tmp_value)])
            @debug "Filtered Results Group:" d1
            @debug "Filtered Ground Truth Group:" d2

            push!(tmp_mse, mean(d1.yhat .- d2.yhat) .^ 2)
            @debug d1, d2

        end
        push!(event_weights, Float64(length(tmp_mse)))
        push!(event_mse, mean(tmp_mse))
    end

    @debug event_mse, event_weights
    # Calculate the mean of the MSE values
    if weighted_mse
        return mean(event_mse, weights(event_weights))
    else
        return event_mse, weights(event_weights)
    end
end