#  # [Cross-validated MANOVA (cvMANOVA)](@id cvmanova_tutorial)
# This tutorial demonstrates how to perform cross-validated MANOVA on a two-level factor.

# # Setup

using Unfold
using UnfoldSim
using UnfoldStats
using StatsModels
using Random
using Statistics
using CairoMakie
using UnfoldMakie
using DataFrames
# Set seed for reproducible cross-validation folds.
Random.seed!(42)

# # Generate two-class EEG data



# Generate simulated multichannel EEG data with two conditions (classes).
sfreq = 100

p1 = (p100(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict())
n1 = (n170(; sfreq = sfreq), @formula(0 ~ 1 + animal), [5, 3], Dict())
p3 = (p300(; sfreq = sfreq), @formula(0 ~ 1 + animal + vegetable), [5, -1, 1], Dict())
design =
    SingleSubjectDesign(;
        conditions = Dict(:animal => ["dog", "cat"], :vegetable => ["tomato", "carrot"]),
        event_order_function = shuffle,
    ) |> x -> RepeatDesign(x, 100)
eeg, evts = UnfoldSim.predef_eeg(
    MersenneTwister(1),
    design,
    LinearModelComponent,
    [p1, n1, p3];
    sfreq,
    return_epoched = true,
    multichannel = true,
    noiselevel = 1.0,
)


fake_times = 1:size(eeg, 2)
# Check the event structure with condition variable
first(evts, 5)

# # Fit an Unfold-Model

# For cvMANOVA we need a overspecified designmatrix. Thus instead of using an intercept, we have one column for each level of any categorical predictor we want to use.
f = @formula 0 ~ 0 + animal + vegetable
contrasts = Dict(
    :animal => StatsModels.FullDummyCoding(),
    :vegetable => StatsModels.FullDummyCoding(),
)

m = fit(UnfoldModel, f, evts, eeg, fake_times; contrasts = contrasts, fit = false)

# As you can see, we have a binary designmatrix now. This is later important to define the contrasts.
modelmatrix(m)[1][1:5, :]

# # Fit with k-fold cross-validation

# We use the cross-validation solver, and importantly, we directly fit the test-set Β estimates as well
cv_solver = Unfold.solver_cv(n_folds = 5, shuffle = true, fit_test = true)
fit!(m, eeg; solver = cv_solver) # note: we could have run the fit(...;solver=...) directly

# # Run cvMANOVA

# Define a two-level contrast to compare the two classes: [-1, 1]
# This means: (class 2) - (class 1)
C = [-1, 1, 0, 0]

# We use only the "baseline" period (first 10 time samples) to estimate a noise covariance.

Y_baseline = eeg[:, 1:10, :]

# Compute cvMANOVA's D for each timepoint 
D_per_fold = cvMANOVA(m, Y_baseline; C = C)

# Aggregate the discriminability statistic across all CV folds.
D_mean = mean(D_per_fold)

# This is the time-resolved discriminability: higher values = better discrimination
let
    f, ax, h = series(reduce(hcat, D_per_fold)', linestyle = :dot)
    lines!(ax, D_mean, color = :black)

    plot_erp!(f[2, 1], subset(coeftable(m), :channel => (x -> x .== 10)))
    f
end

# ## Cross Decoding
# Next we will check how well we can decode animal based on vegetable, and vice versa
C_animal = [-1, 1, 0, 0]
C_veggie = [0, 0, -1, 1]
D_animal = mean(cvMANOVA(m, Y_baseline; C = C_animal))
D_veggie = mean(cvMANOVA(m, Y_baseline; C = C_veggie))
D_cross = mean(cvMANOVA(m, Y_baseline; C = C_animal, C_test = C_veggie))

let
    f = Figure()
    ax = Axis(f[1, 1])
    h_a = lines!(ax, D_animal, label = "animal")
    h_v = lines!(ax, D_veggie, label = "veggie")
    h_c = lines!(ax, D_cross, label = "cross")
    Legend(f[1, 2], ax)
    f
end

# Cross decoding is only possible where veggie and animal share some representation, e.g. at the "p300" time window, but not at the "n170" time window, which was simulated specific to animals.

## Time Generalization
import LinearAlgebra: diag

D_temporal = mean(
    cvMANOVA(
        m,
        Y_baseline;
        C = C_animal,
        C_test = C_veggie,
        temporal_generalization = true,
    ),
)
let
    f, ax, h = heatmap(D_temporal, axis = (; aspect = DataAspect()))
    lines!(ax, [0, size(D_temporal, 1)], [0, size(D_temporal, 2)]) # diag
    lines(f[1, 2], diag(D_temporal), label = "diag(temp-gen)")
    lines!(D_cross, linestyle = :dash, label = "vector")
    Legend(f[1, 3], current_axis())
    Label(f[1, 1, TopLeft()], "A)")
    Label(f[1, 1, TopLeft()], "B)")
    f
end
# A) shows the temporal generalization matrix, where the x-axis is the training time and the y-axis is the testing time. Here one can see that training at sample 20 allows decoding at sample 30
# B) The diagonal of this matrix is equivalent to the "normal" way as discussed before