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
nothing##hide
# # Generate two-class EEG data

# Generate simulated multichannel EEG data with two conditions (classes).
eeg, evts =
    UnfoldSim.predef_eeg(; return_epoched = true, multichannel = true, noiselevel = 1.0)
fake_times = 1:size(eeg, 2)
nothing ##hide
# Check the event structure with condition variable
first(evts, 5)

# # Fit an Unfold-Model

# For cvMANOVA we need a overspecified designmatrix. Thus instead of using an intercept, we have one column for each level of any categorical predictor we want to use.
f = @formula 0 ~ 0 + condition
contrasts = Dict(:condition => StatsModels.FullDummyCoding())

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
C = [-1, 1]

# We use only the "baseline" period (first 10 time samples) to estimate a noise covariance.

Y_baseline = eeg[:, 1:10, :]
nothing##hide
# Compute cvMANOVA's D for each timepoint 
D_per_fold = cvMANOVA(m, Y_baseline; C = C)

# Aggregate the discriminability statistic across all CV folds.
D_mean = mean(D_per_fold)
nothing ##hide
# This is the time-resolved discriminability: higher values = better discrimination
f, ax, h = series(reduce(hcat, D_per_fold)', linestyle = :dot)
lines!(ax, D_mean, color = :black)

plot_erp!(f[2, 1], subset(coeftable(m), :channel => (x -> x .== 1)))
f
