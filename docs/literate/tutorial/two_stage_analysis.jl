# # Two-stage EEG analysis using Unfold & UnfoldStats

# ## 0. Setup
# Load required packages. 
using UnfoldStats
using Unfold
using DataFrames
using Chain
using Statistics
using HypothesisTests
using CairoMakie
using UnfoldMakie

# ## 1. Simulate data
# This section can be skipped, if one already has (real) data that they want to analyse.

# ```@raw html
# <details>
# <summary>Click to expand the simulation details</summary>
# ```

using UnfoldSim
using StableRNGs

design = MultiSubjectDesign(
    n_subjects = 20,
    n_items = 100,
    items_between = Dict(:continuous => range(0, 2, length = 10)),
)

β = [1, 1, 0.5, 0.2]
σs = Dict(:subject => [1, 0.1, 0.1, 0.1])

signal = MixedModelComponent(;
    basis = UnfoldSim.hanning(50),
    formula = @formula(
        0 ~
            1 +
            continuous +
            continuous^2 +
            continuous^3 +
            (1 + continuous + continuous^2 + continuous^3 | subject)
    ),
    β = β,
    σs = σs,
)

hart = headmodel(type = "hartmut")
signal_multichannel = MultichannelComponent(signal, hart => "Left Postcentral Gyrus")

onset = UniformOnset(; width = 50, offset = 60)
noise = PinkNoise(; noiselevel = 5)

# ```@raw html
# </details >
# ```

data, events = simulate(
    StableRNG(1),
    design,
    signal_multichannel,
    onset,
    noise,
    return_epoched = false,
);

# ```@raw html
# <details>
# <summary>Click to expand event data frame </summary>
# ```
first(events, 12)
# ```@raw html
# </details >
# ```

# ## 2. Fit an Unfold model for each subject (first stage)
# In the first stage, we fit an Unfold model for each subject separately.

## Specify a temporal basis function
basisfunction = firbasis((-0.1, 0.7), 100)

## Specify the model formula
formula = @formula 0 ~ 1 + spl(continuous, 4)

## Combine basisfunction and formula in an event dict
event_dict = Dict(Any => (formula, basisfunction))

subject_list = unique(events.subject)
model_list = UnfoldLinearModelContinuousTime[]

## Slice the data by its last dimension (i.e. the subject dimension)
data_slices = eachslice(data, dims = ndims(data))

for s = 1:size(data, ndims(data))
    m = fit(
        UnfoldModel,
        event_dict,
        subset(events, :subject => ByRow(==(subject_list[s]))),
        data_slices[s],
    )
    push!(model_list, m)
end

models = DataFrame(subject = subject_list, unfoldmodel = model_list);
size(models)

# As a result, we get a dataframe which contains an Unfold model for each subject.

# ## 3. Compute and visualize the marginal effects

# In the next step, we will compute the marginal effects of the `continuous` predictor i.e. how the prediction changes
# for different levels of this predictor. First, we compute the marginal effects separately for each subject.
# Then we aggregate them over subjects using the `mean`.

## Specify the predictors of interest
predictor_dict = Dict(:continuous => range(0, 2, length = 3))

## Compute the marginal effects for all subjects separately
effects_all_subjects = combine(
    groupby(models, :subject),
    :unfoldmodel => (m -> effects(predictor_dict, m[1])) => AsTable,
)

## Aggregate the marginal effects (per event type, time point and channel) over subjects
## using the mean as aggregation function
aggregated_effects = @chain effects_all_subjects begin
    groupby([:basisname, :channel, collect(keys(predictor_dict))..., :time])
    combine(:yhat .=> [x -> mean(skipmissing(x))] .=> Symbol("yhat_", mean))
end;
first(aggregated_effects, 5)

# Next, we want to visualize the marginal effects in a topoplot series. For this we first need to extract the electrode
# positions for our simulated data (or your real data) and project them from 3d to 2d. 

## Extract the electrode positions for the simulated data from the headmodel
## and project them from 3d to 2d.
pos3d = hart.electrodes["pos"];
pos2d = to_positions(pos3d')
pos2d = [Point2f(p[1] + 0.5, p[2] + 0.5) for p in pos2d];

# For visualization, we use the [`plot_topoplotseries` function](https://unfoldtoolbox.github.io/UnfoldMakie.jl/dev/generated/tutorials/topoplotseries/) from `UnfoldMakie`.
# In the topoplot series, the time in seconds (binned in time windows) is represented from left to right. 
# The rows (from top to bottom) represent the marginal effects for different levels of the predictor `continuous`.

## Set the size of the time bins for the topoplot series
bin_size = 0.1

f_effects = Figure(size = (1200, 600))
tp_effects = plot_topoplotseries!(
    f_effects,
    aggregated_effects,
    bin_size,
    positions = pos2d,
    mapping = (; y = :yhat_mean, row = :continuous),
    visual = (; enlarge = 0.6, label_scatter = false, colorrange = (-3, 3)),
)

ax = current_axis()
linkaxes!(tp_effects.content[1:end-2]...)
xlims!(ax, 0, 0.9)
ylims!(ax, 0, 0.9)
current_figure()

# In the time windows `[0.1, 0.2)`, `[0.2, 0.3)` and `[0.3, 0.4)`, one can see the effect of `continuous` that we simulated.
# In the next section, we want to quantify and test this effect.

# ## 4. Extract coefficients & conduct Hotelling's T² tests

# We will use a one-sample Hotelling's T² test to test whether at least one of the spline coefficients is different from 0.
# First, we extract the spline coefficients (for the `continuous` predictor variable) for all subjects. Then we conduct
# a Hotelling's T² test separately for each channel and time point and extract the correspoding p-value. 

## Extract the spline coefficients for the `continuous` predictor variable
## from the Unfold models of all subjects
basisname = unique(effects_all_subjects.basisname)
coefs = extract_coefs(models.unfoldmodel, :continuous, basisname)

## Conduct a one-sample Hotelling's T² test separately for all channels and time points
## and compute a p-value. We compare the spline coefficients vector against 0.
p_values =
    mapslices(c -> pvalue(OneSampleHotellingT2Test(c', [0, 0, 0])), coefs, dims = (3, 4)) |>
    x -> dropdims(x, dims = (3, 4));

times = unique(effects_all_subjects.time)
channels = axes(p_values, 1)

## For visualization purposes, save the p_values in a data frame together with time and channel
p_values_df = DataFrame(
    channel = repeat(channels, outer = length(times)),
    time = repeat(times, inner = length(channels)),
    p_values = p_values[:],
);

first(p_values_df, 5)

# As a last step, we visualize the p-values in a topoplot series.
bin_size = 0.1

f_pvalues = Figure(size = (1200, 200))
tp_pvalues = plot_topoplotseries!(
    f_pvalues,
    p_values_df,
    bin_size,
    positions = pos2d,
    mapping = (; y = :p_values),
    visual = (; enlarge = 0.6, label_scatter = false, colorrange = (0, 0.1)),
    colorbar = (; label = "p-value"),
)

ax = current_axis()
linkaxes!(tp_pvalues.content[1:end-2]...)
xlims!(ax, 0, 0.9)
ylims!(ax, 0, 0.9)
current_figure()

# !!! note
#     To decide whether to consider the effect statistically significant and to correct for multiple comparisons
#     (due to different time points and channels), one could conduct a cluster-permutation test using the Hotelling's T² 
#     values as the test statistic.
