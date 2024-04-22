module UnfoldStatsMixedModelsExt
using Unfold
using UnfoldStats
using MixedModels


lmm_ext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)

if isnothing(lmm_ext)
    error("Something went wrong with getting the Unfold UnfoldMixedModelsExt extension")
end
# Currently, `extract_coefs` is not implemented for mixed-effects models
UnfoldStats.extract_coefs(
    model::Union{
        lmm_ext.UnfoldLinearMixedModel,
        lmm_ext.UnfoldLinearMixedModelContinuousTime,
    },
    predictor,
    basisname,
) = throw(
    "The `extract_coefs` function is currently not implemented for mixed-effects models.",
)
end