# !!! important
#       This functionality is not extensively tested. While it returns reasonable results, we haven't written any unit-tests, nor tested the type-1 error probability yet


# ## 0. Setup
# If you want to do LMM + cluster permutations, you at least need these packages (more for the simulation below)
using UnfoldMixedModels
using MixedModelsPermutations, ClusterDepth
using UnfoldStats

# ## 1. Simulate data
# This section can be skipped, if one already has (real) data that they want to analyse.



# ```@raw html
# <details>
# <summary>Click to expand the simulation details</summary>
# ```
using UnfoldSim
using StatsModels
using Random

srate = 25
design = MultiSubjectDesign(;
    n_subjects = 30,
    n_items = 40,
    items_between = Dict(:stimtype => ["car", "face"]),
)
#both_within = Dict(:condition=>["scrambled","intact"]))
contrasts = Dict(:stimtype => DummyCoding())
p1 = MixedModelComponent(;
    basis = UnfoldSim.p100(; sfreq = srate),
    formula = @formula(dv ~ 1 + (1 | subject) + (1 | item)),
    β = [5.0],
    σs = Dict(:subject => [0.0], :item => [0.0]),
    contrasts = contrasts,
);

n1 = MixedModelComponent(;
    basis = UnfoldSim.n170(; sfreq = srate),
    formula = @formula(dv ~ 1 + stimtype + (1 + stimtype | subject) + (1 | item)),
    β = [1.0, 4], # n170-basis is negative
    σs = Dict(:subject => [2.0, 0.25], :item => [0.25]),
    contrasts = contrasts,
);

p3 = MixedModelComponent(;
    basis = UnfoldSim.p300(; sfreq = srate),
    formula = @formula(dv ~ 1 + (1 | subject) + (1 + stimtype | item)),
    β = [4.0],
    σs = Dict(:subject => [1.0], :item => [0.5, 2]),
    contrasts = contrasts,
);



data_e, events = UnfoldSim.simulate(
    design,
    [p1, n1, p3],
    UniformOnset(srate * 2, 10),
    PinkNoise(; noiselevel = 1);
    return_epoched = true,
)
times = range(-0.1, 0.5, length = size(data_e, 1))
data_e = reshape(data_e, 1, size(data_e, 1), :)
nothing ##hide
# ```@raw html
# </details >
# ```

# ## 2. Fit mass-univariate LMMs
# We have some typical experimental data with subject and item effects. Item refer to stimuli here, based on our `stimtype` condition these are either different `cars` or `faces`.
m = fit(
    UnfoldModel,
    [
        Any => (
            @formula(0 ~ 1 + stimtype + (1 + stimtype | item) + (1 + stimtype | subject)),
            times,
        ),
    ],
    events,
    data_e,
);


# ## 3. Cluster permutation test
# If we would run a statistical test on each time-point separately, we would greatly inflate the type-1 error, reaching significance on any each sample much higher than the assumed α=0.05.
# One solution are cluster permutation test, where we instead test for clustersizes of connected significant clusters. In "classical" two-stage testing, such a permutation test is straight forward. But for LMMs we have to think of something more clever, as it is not directly clear how to permute if both subject and item effects exist (you gonna break the relation between the two). We did that in `MixedModelsPermutations` and can apply this strategy to EEG data as well.`

# select the fixed-effects coefficient to test (`stimtype`)
coefficient = 2;

# call the permutation test
# !!! note
#      This interface is very likely to change in the future
pvalue(
    MersenneTwister(1),
    m,
    data_e,
    coefficient;
    n_permutations = 20,
    clusterforming_threshold = 1.8,
)

