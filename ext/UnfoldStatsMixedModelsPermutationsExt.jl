module UnfoldStatsMixedModelsPermutationsExt

using UnfoldMixedModels
using UnfoldStats
#using MixedModels
using MixedModelsPermutations
using ClusterDepth
using Logging
using Random



"""
        UnfoldStats.pvalue(rng,model::UnfoldLinearMixedModel,data,coefficient;kwargs)

Calculate pvalues of an UnfoldLinearMixedModel of various types, and various multiple-comparison corrections.

While We write various here, we mean that currently only `ClusterDepth` is actually implemented.

# Arguments
- `rng::AbstractRNG`: An RNG-generator for reproducibility
- `model::UnfoldLinearMixedModel`: The fitted UnfoldLinearMixedModel to calculate the permutations of
- `data::AbstractArray{<:Real,3}`: The data with ch x time x trials
- `coefficient:Int`: The coefficient to test the null of


# Keyword arguments
- `type::String = "clusterdepth": Specify the type. Only "clustedepth" is currently implemented. "FDR", or "Walds'z uncorrected" will follow
- `clusterforming_threshold::Float`: The threshold used to calculate what is a cluster and what not. We use the `abs` around observed & permuted data.
- `n_permutations::Int = 500`: Number of permutations. Based on other permutation work, 500 should be a reasonable number. Methods based on e.g. tail-approximation are not yet available
- `lmm_statistic = :z`: What statistic to extract, could also be e.g. `β`
- `time_selection = 1:size(data, 2)`: Possibility to calculate permutations on a subset of time points


# Returns
- `result::Array` : An array of ch x time pvalues
"""
function UnfoldStats.pvalue(
    rng,
    model::UnfoldLinearMixedModel,
    data,
    coefficient;
    type = "clusterdepth",
    clusterforming_threshold,
    kwargs...,
)
    if type != "clusterdepth"
        error("other types (e.g. FDR currently not implemented")
    elseif type == "clusterdepth"
        pvals = lmm_clusterdepth(
            rng,
            model,
            data,
            coefficient;
            clusterforming_threshold,
            kwargs...,
        )
        return reshape(pvals, size(data)[1:end-1]...)
    end
end



"""
        lmm_clusterdepth(rng,model,data,coefficient;kwargs...)

Entry point to calculate cluster-corrected pvalues for LMMs. Extracts the observed values, calculates permutations, and finally returns Troendle-corrected pvalues.


# Arguments (if needed)
- `rng::AbstractRNG`: An RNG-generator for reproducibility
- `model::UnfoldLinearMixedModel`: The fitted UnfoldLinearMixedModel to calculate the permutations of
- `data::AbstractArray{<:Real,3}`: The data with ch x time x trials
- `coefficient:Int`: The coefficient to test the null of

# Keyword arguments
See `lmm_clusterdepth_pvalues`, `lmm_permutations` and `get_lmm_statistic`.

# Returns
- `result::Vector` : cluster corrected pvalues

"""
function lmm_clusterdepth(
    rng,
    model,
    data,
    coefficient;
    lmm_statistic = :z,
    clusterforming_threshold,
    kwargs...,
)
    permuted = lmm_permutations(rng, model, data, coefficient; kwargs...)
    observed = get_lmm_statistic(model, coefficient, lmm_statistic)
    return lmm_clusterdepth_pvalues(rng, observed, permuted; clusterforming_threshold)
end




"""
        lmm_clusterdepth_pvalues(rng, observed, permuted; clusterforming_threshold)

Returns per-sample two-sided pvalues based on observed data and permuted, after applying a clusterforming_threshold

# Arguments
- `rng::AbstractRNG`: An RNG-generator for reproducibility
- `observed::Vector`: The observed-data, linearized (that is, instead of ch x time, we do need the `x[:]` equivalent)
- `permuted::Array`: The ch x time x n_permutations matrix under the H₀

# Keyword arguments
- `clusterforming_threshold::Float`: The threshold used to calculate what is a cluster and what not. We use the `abs` around observed & permuted data.

# Returns
- `result::Vector` : pvalues in the same size as observed
"""
function lmm_clusterdepth_pvalues(rng, observed, permuted; clusterforming_threshold)

    # we need global variables here (yes, sorry...), because instead of actually
    # letting ClusterDepth do the permutation, we just have to index the already
    # permuted data given in the function (`permuted`)
    global n_permutation_count
    n_permutation_count = 0
    function _fake_permutation_fun(r, data)
        global n_permutation_count
        n_permutation_count = n_permutation_count + 1
        return permuted[:, :, n_permutation_count]
    end
    J_tuple = ClusterDepth.perm_clusterdepths_both(
        rng,
        abs.(permuted),
        _fake_permutation_fun,
        clusterforming_threshold;
        statfun = x -> abs.(x),
        nₚ = size(permuted, 3),
        sidefun = x -> x
    )

    pvals = ClusterDepth.pvals(abs.(observed), J_tuple, clusterforming_threshold)

end


"""
    lmm_permutations(rng::AbstractRNG,model::UnfoldLinearMixedModel,data::AbstractArray{<:Real,3},coefficient::Int;kwargs...)

Calculates permutations of a UnfoldLinearMixedModel object

# Arguments
- `rng::AbstractRNG`: An RNG-generator for reproducibility
- `model::UnfoldLinearMixedModel`: The fitted UnfoldLinearMixedModel to calculate the permutations of
- `data::AbstractArray{<:Real,3}`: The data with ch x time x trials
- `coefficient:Int`: The coefficient to test the null of

# Keyword arguments
- `n_permutations::Int = 500`: Number of permutations. Based on other permutation work, 500 should be a reasonable number. Methods based on e.g. tail-approximation are not yet available
- `lmm_statistic = :z`: What statistic to extract, could also be e.g. `β`
- `time_selection = 1:size(data, 2)`: Possibility to calculate permutations on a subset of time points

# Returns
- `result::Array` : Returns the full permutation Array with ch x timepoints x permutations for the coefficient-of-interest
"""
function lmm_permutations(
    rng::AbstractRNG,
    model,
    data::AbstractArray{<:Real,3},
    coefficient::Int;
    n_permutations = 500,
    lmm_statistic = :z,
    time_selection = 1:size(data, 2),
)
    permdata = Array{Float64}(undef, size(data, 1), length(time_selection), n_permutations)


    Xs = UnfoldMixedModels.prepare_modelmatrix(model)

    mm_outer = UnfoldMixedModels.LinearMixedModel_wrapper(
        Unfold.formulas(model),
        data[1, 1, :],
        Xs,
    )
    mm_outer.optsum.maxtime = 0.1 # 

    chIx = 1 # for now we only support 1 channel anyway
    #
    #p = Progress(length(time_selection))

    @assert(
        coefficient <= size(coef(mm_outer))[end],
        "chosen coefficient was larger than available coefficients"
    )
    #Threads.@threads for tIx =1:length(time_selection)
    #@showprogress "Processing Timepoints" for 
    for chIx = 1:size(data, 1)
        #Threads.@threads 
        for tIx = 1:length(time_selection)

            # splice in the correct dataa for residual calculation
            mm = deepcopy(mm_outer)
            mm.y .= data[chIx, time_selection[tIx], :]

            # set the previous calculated model-fit
            θ = Vector(modelfit(model).θ[time_selection[tIx]])
            @debug size(θ)
            MixedModels.updateL!(MixedModels.setθ!(mm, θ))

            # get the coefficient 
            H0 = coef(mm)
            # set the one of interest to 0
            H0[coefficient] = 0
            # run the permutation

            permutations = undef
            Logging.with_logger(NullLogger()) do   # remove NLopt warnings  
                permutations = permutation(
                    deepcopy(rng), # important here is to set the same seed to keep flip all time-points the same
                    n_permutations,
                    mm;
                    β = H0,
                    progress = false,
                    #blup_method = MixedModelsPermutations.olsranef,
                ) # constant rng to keep autocorr & olsranef for singular models
            end

            # extract the test-statistic

            permdata[chIx, tIx, :] =
                get_lmm_statistic(model, permutations, coefficient, lmm_statistic)

            #next!(p)
        end # end for
    end
    return permdata
end



"""
    get_lmm_statistic(model::UnfoldLinearMixedModel, coefficient::Int, lmm_statistic)
    get_lmm_statistic(model,permutations::MixedModelFitCollection, coefficient, lmm_statistic)
    
Returns the field `lmm_statistic` of the `coefpvalues` table output of the `permutations`-FitCollection only of the coefficient `coefnames(formulas(model))[1][coefficient]`

# Arguments
- `model`: typically an `UnfoldMixedModel` or similar
- `permutations::MixedModelFitCollection`: the output of `modelfit(model)`
- `coefficient::Int`: the coefficient to choose (a fixed effect)
- `lmm_statistic::Symbol`: The statistic to extact, tyically either `β` or `z`


# Returns
- `result::Vector` : A vector of the extracted `lmm_statistic`

"""
function get_lmm_statistic(
    model,
    permutations::MixedModelsPermutations.MixedModels.MixedModelFitCollection,
    coefficient::Int,
    lmm_statistic,
)
    [
        getproperty(m, lmm_statistic) for m in permutations.coefpvalues if
        String(m.coefname) == Unfold.coefnames(Unfold.formulas(model))[1][coefficient]
    ]

end
function get_lmm_statistic(model::UnfoldLinearMixedModel, coefficient::Int, lmm_statistic)
    return get_lmm_statistic(model, modelfit(model), coefficient, lmm_statistic)
    #    r = coeftable(m)
    #    r = subset(r, :group => (x -> isnothing.(x)), :coefname => (x -> x .!== "(Intercept)"))
    #    tvals = abs.(r.estimate ./ r.stderror)
    #    return tvals
end
end
