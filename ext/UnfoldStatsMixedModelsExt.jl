module UnfoldStatsMixedModelsExt
using Unfold
using UnfoldStats

lmm_ext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)
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