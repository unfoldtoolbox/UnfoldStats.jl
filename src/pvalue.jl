# needed for pvalues in the extension
#@deprecate pvalues(args...) pvalue(args...)
StatsAPI.pvalue(args...; kwargs...) =
    error("not implemented, see `methods(pvalues)` for a list")