module UnfoldStatsMixedModelsPermutationsExt
using Unfold
import Unfold: pvalues
using UnfoldStats
#using MixedModels
using MixedModelsPermutations
using ClusterDepth
using Logging
using Random
const MixedModels = MixedModelsPermutations.MixedModels
LMMext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)


function Unfold.pvalues(
    rng,
    model::LMMext.UnfoldLinearMixedModel,
    data,
    coefficient;
    type = "clusterdepth",
    clusterforming_threshold,
    kwargs...,
)
    if type != "clusterdepth"
        error("other types (e.g. FDR currently not implemented")
    elseif type == "clusterdepth"
        return lmm_clusterdepth(
            rng,
            model,
            data,
            coefficient;
            clusterforming_threshold,
            kwargs...,
        )
    end
end




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
    return lmm_clusterdepth_pvalues(
        rng,
        observed,
        permuted;
        clusterforming_threshold,
        kwargs...,
    )
end
function lmm_clusterdepth_pvalues(
    rng,
    observed,
    permuted;
    clusterforming_threshold,
    kwargs...,
)

    # we need global variables here (yes, sorry...), because instead of actually
    # letting ClusterDepth do the permutation, we just have to index the already
    # permuted data given in the function (`permuted`)
    global n_permutation_count
    n_permutation_count = 0
    function _fake_permutation_fun(r, data)
        global n_permutation_count
        n_permutation_count = n_permutation_count + 1
        return permuted[:, n_permutation_count]
    end
    J_tuple = ClusterDepth.perm_clusterdepths_both(
        rng,
        abs.(permuted),
        _fake_permutation_fun,
        clusterforming_threshold;
        statfun = x -> abs.(x),
        nₚ = size(permuted, 2),
    )

    pvals = ClusterDepth.pvals(abs.(observed), J_tuple, clusterforming_threshold)

end

function lmm_permutations(
    rng::AbstractRNG,
    model,
    data::AbstractArray{<:Real,3},
    coefficient::Int;
    n_permutations = 500,
    lmm_statistic = :z,
    time_selection = 1:size(data, 2),
)
    permdata = Matrix{Float64}(undef, length(time_selection), n_permutations)


    Xs = LMMext.prepare_modelmatrix(model)

    mm_outer = LMMext.LinearMixedModel_wrapper(Unfold.formulas(model), data[1, 1, :], Xs)
    mm_outer.optsum.maxtime = 0.1 # 

    chIx = 1 # for now we only support 1 channel anyway
    #
    #p = Progress(length(time_selection))
    #Threads.@threads for tIx =1:length(time_selection)
    #@showprogress "Processing Timepoints" for 
    Threads.@threads for tIx = 1:length(time_selection)

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
                hide_progress = true,
                #blup_method = MixedModelsPermutations.olsranef,
            ) # constant rng to keep autocorr & olsranef for singular models
        end

        # extract the test-statistic

        permdata[tIx, :] =
            get_lmm_statistic(model, permutations, coefficient, lmm_statistic)

        #next!(p)
    end # end for
    return permdata
end

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
function get_lmm_statistic(
    model::LMMext.UnfoldLinearMixedModel,
    coefficient::Int,
    lmm_statistic,
)
    return get_lmm_statistic(model, modelfit(model), coefficient, lmm_statistic)
    #    r = coeftable(m)
    #    r = subset(r, :group => (x -> isnothing.(x)), :coefname => (x -> x .!== "(Intercept)"))
    #    tvals = abs.(r.estimate ./ r.stderror)
    #    return tvals
end
end