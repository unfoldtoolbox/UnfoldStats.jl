# Helper functions
const BSplineTerm = Base.get_extension(Unfold, :UnfoldBSplineKitExt).BSplineTerm

extract_symbol(t::AbstractTerm) = t.sym

extract_symbol(t::InterceptTerm) = "(Intercept)"
extract_symbol(t::BSplineTerm) = repeat([t.term.sym], t.df - 1)
extract_symbol(t::InteractionTerm) = extract_symbol.(t.terms)
extract_symbol(t::FunctionTerm) = extract_symbol.(t.args)

extract_symbol(t::MatrixTerm) = extract_symbol(t.terms)
extract_symbol(t::Unfold.TimeExpandedTerm) =
    repeat(extract_symbol(t.term), inner = length(Unfold.colnames(t.basisfunction)))
extract_symbol(f::FormulaTerm) = vcat(extract_symbol.(f.rhs)...)
extract_symbol(t::Vector) = vcat(extract_symbol.(t)...)
extract_symbol(t::Tuple) = vcat(extract_symbol.(t)...)


contained_or_equal(p, e) = (p == e)
contained_or_equal(p, e::Tuple) = (p in e)

get_predictor_string(p::Symbol) = ":$p"
get_predictor_string(p::String) = "\"$p\""
get_predictor_string(p::Tuple) = "$p"

function get_basisnames(model::UnfoldLinearModel)
    # Extract the event names from the design
    design_keys = keys(Unfold.design(model))

    # Create a list of the basis names corresponding to each model term
    basisnames = String[]
    for (ix, event) in enumerate(design_keys)
        push!(basisnames, repeat(["event: $(event)"], size(modelmatrix(model)[ix], 2))...)
    end
    return basisnames
end

get_basisnames(model::UnfoldLinearModelContinuousTime) =
    Unfold.extract_coef_info(Unfold.get_coefnames(model), 1)


function extract_coefs(model::UnfoldModel, predictor, basisname)

    # Get vector with underlying predictor variable (symbol) for all coefficients
    symbols = extract_symbol(Unfold.formula(model))

    # Check whether `predictor` is a predictor in the model
    if predictor ∉ symbols
        # TODO: Interactions will be listed separately at the moment, maybe don't list them
        allowed_predictors = join(get_predictor_string.(unique(symbols)), ", ")

        throw(
            ArgumentError(
                "The given predictor $(get_predictor_string(predictor)) was not found in the model. Possible predictors are: $allowed_predictors.",
            ),
        )
    end

    basisname_list = get_basisnames(model)

    # Check whether given `basisname` is in the basisname list of the model
    if basisname ∉ basisname_list
        allowed_basisnames = join(["\"$b\"" for b in unique(basisname_list)], ", ")

        throw(
            ArgumentError(
                "The given basisname \"$basisname\" was not found in the model. Possible basisnames are: $allowed_basisnames.",
            ),
        )
    end

    # Create a boolean mask which is true for the model coefficients that belong to the given predictor
    mask_predictor = contained_or_equal.(predictor, symbols)

    # Create a boolean mask which is true for the model coefficients that belong to the given basis name
    mask_basisfunction = basisname_list .== basisname

    mask_combined = mask_predictor .* mask_basisfunction

    # Check whether the given combination between predictor variable and basisname exist in the model
    if sum(mask_combined) == 0
        # TODO: Is `ArgumentError` the right exception to use?
        throw(
            ArgumentError(
                "The given predictor $(get_predictor_string(predictor)) does not exist for the given basisname \"$basisname\".",
            ),
        )
    end

    # Extract the requested coefficients from the coefficient array
    if typeof(model) == UnfoldLinearModel
        coef_subset = coef(model)[:, :, mask_combined]

    elseif typeof(model) == UnfoldLinearModelContinuousTime

        n_coefs = length(
            unique(Unfold.extract_coef_info(Unfold.get_coefnames(model), 2)[mask_combined]),
        )

        coef_subset_temp = @view(coef(model)[:, mask_combined])
        # Reshape the coefficient array such that time points and coefficient types get separate dimensions
        coef_subset = reshape(coef_subset_temp, size(coef_subset_temp, 1), :, n_coefs)
    else
        throw("Not implemented.")
    end

    return coef_subset#::Array{<:Union{<:Missing,<:Float64},3}
end

# Currently, `extract_coefs` is not implemented for mixed-effects models
extract_coefs(
    model::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
    predictor,
    basisname,
) = throw(
    "The `extract_coefs` function is currently not implemented for mixed-effects models.",
)

function extract_coefs(models::Vector{<:UnfoldModel}, predictor, basisname)

    # Extract the coefficients for all subjects
    coefs_vector = extract_coefs.(models, predictor, basisname)

    # Check that all coefficient arrays have the same size
    @assert length(unique(size.(coefs_vector))) == 1

    # Concatenate the coefficients to one array
    # Dimensions: (channels x times x coefficients x subjects)
    coefs_all_subjects = cat(coefs_vector..., dims = ndims(coefs_vector[1]) + 1)#::Array{<:Union{<:Missing,<:Float64},4}

    return coefs_all_subjects
end