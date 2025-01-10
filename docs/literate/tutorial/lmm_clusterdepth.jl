
# get some data
using UnfoldSim
using Unfold
using MixedModelsPermutations, ClusterDepth # both necessary to activate correct extension!
using UnfoldStats
using StatsModels
using Random
using MixedModels

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
) # 18
times = range(-0.1, 0.5, length = size(data_e, 1))
data_e = reshape(data_e, 1, size(data_e, 1), :)
#events.latency .+= repeat(range(0,length=size(data,2),step=size(data,1)),inner=size(events[events.subject.=="S01",:],1))



# # Fit LMM 
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


# # Cluster Permute :)
coefficient = 2
pvalues(
    MersenneTwister(1),
    m,
    data_e,
    coefficient;
    n_permutations = 10,
    clusterforming_threshold = 1.8,
)