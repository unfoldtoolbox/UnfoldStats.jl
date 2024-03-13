var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = UnfoldStats","category":"page"},{"location":"#UnfoldStats","page":"Home","title":"UnfoldStats","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for UnfoldStats.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [UnfoldStats]","category":"page"},{"location":"#UnfoldStats.contained_or_equal-Tuple{Any, Any}","page":"Home","title":"UnfoldStats.contained_or_equal","text":"contained_or_equal(p, e)\n\nTest if p equals e or whether e contains p if e is a tuple.\n\n\n\n\n\n","category":"method"},{"location":"#UnfoldStats.extract_coefs-Tuple{Unfold.UnfoldModel, Any, Any}","page":"Home","title":"UnfoldStats.extract_coefs","text":"extract_coefs(model::UnfoldModel, predictor, basisname)\n\nReturn the coefficients of an Unfold model for a certain predictor and basisname.\n\nFor extracting the terms of a predictor variable predictor must be a symbol e.g. :continuous. For extracting the intercept predictor should be a String, i.e. \"(Intercept)\".\n\nbasisname must match one of the basis names which can be found in coeftable(model).\n\nNote: If a predictor variable has more than one term in the formula (e.g. a spline set, a categorical variable with several levels or an interaction), the coefficients for all terms are returned.\n\nThe dimensions of the returned coefficients are channel x times x coefficients.\n\n\n\n\n\n","category":"method"},{"location":"#UnfoldStats.extract_coefs-Tuple{Vector{<:Unfold.UnfoldModel}, Any, Any}","page":"Home","title":"UnfoldStats.extract_coefs","text":"extract_coefs(models::Vector{<:UnfoldModel}, predictor, basisname)\n\nWhen applied to a vector of Unfold models, extracts the coefficients (matching the predictor and basisname) for all models (usually subjects) and concatenates them.\n\nThe dimensions of the returned coefficients are channel x times x coefficients x subjects.\n\n\n\n\n\n","category":"method"},{"location":"#UnfoldStats.extract_symbol-Tuple{StatsModels.AbstractTerm}","page":"Home","title":"UnfoldStats.extract_symbol","text":"extract_symbol(t::AbstractTerm)\n\nReturn the symbol(s) underlying a term from a model formula, repeated by their actual coefficient number (after StatsModels.apply_schema).\n\nExamples\n\njulia> f = @formula 0 ~ 1 + spl(continuous, 4) + continuous + condition + pet + condition & pet\njulia> ... # apply schema using an event dataframe, according to StatsModels\njulia> extract_symbol(f)\n8-element Vector{Any}:\n \"(Intercept)\"\n :continuous\n :continuous\n :continuous\n :continuous\n :condition\n :pet\n (:condition, :pet)\n\nWe get the actual symbols of each predictor - this is different to a function that would return the symbol for each term, which would be [\"(Intercept)\", :continuous,:continuous,:condition,:pet,(:condition,:pet) ]\n\nThe difference between those two cases would get even more stark, if a basisfunction is in play as it timeexpand terms into many more predictors.\n\n\n\n\n\n","category":"method"},{"location":"#UnfoldStats.get_basisnames-Tuple{Unfold.UnfoldLinearModel}","page":"Home","title":"UnfoldStats.get_basisnames","text":"get_basisnames(model::UnfoldModel)\n\nReturn the basisnames for all predictor terms as a vector.\n\nThe returned vector contains the name of the event type/basis, repeated by their actual coefficient number (after StatsModels.apply_schema). If a model has more than one event type (e.g. stimulus and fixation), the vectors are concatenated.\n\n\n\n\n\n","category":"method"},{"location":"#UnfoldStats.get_predictor_string-Tuple{Symbol}","page":"Home","title":"UnfoldStats.get_predictor_string","text":"get_predictor_string(p)\n\nReturn string representation based on the type of p.\n\nThis function is used for a useful display of variables e.g. in an error message.\n\nExamples\n\njulia> UnfoldStats.get_predictor_string(:condition)\n\":condition\"\n\n\n\n\n\n","category":"method"}]
}
