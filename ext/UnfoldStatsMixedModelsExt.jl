module UnfoldStatsMixedModelsExt
using UnfoldMixedModels
using UnfoldStats


# Currently, `extract_coefs` is not implemented for mixed-effects models
UnfoldStats.extract_coefs(
    model::Union{
        UnfoldLinearMixedModel,
        UnfoldLinearMixedModelContinuousTime,
    },
    predictor,
    basisname,
) = throw(
    "The `extract_coefs` function is currently not implemented for mixed-effects models.",
)
    
end