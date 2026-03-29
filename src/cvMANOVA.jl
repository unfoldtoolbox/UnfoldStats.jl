
"""
    calculate_Σinv(Ξ;λ=0.)
Calculate the regularized inverse of the noise covariance matrix from the residuals Ξ, with optional regularization parameter λ.
Regularization towards the diagonal is used, the Ledoit-Wolf shrinkage method.

# Arguments
- `Ξ`: The residuals, used to estimate the noise covariance (ch x time)
- `λ`: Regularization parameter for the noise covariance estimation (default: 0.0 - recommended to be quite small in 2014 paper)
"""
function calculate_Σinv(Ξ; λ = 0.0)
    Σ = Ξ'Ξ # error covariance matrix

    Σ_reg = (1 - λ) * Σ + λ * Diagonal(diag(Σ)) # regularized error covariance matrix
    Σinv = Σ_reg^-1  # invert
end

function cvmanova_D(CCXXCC, Σinv, B_train, B_test; scaling_factor = 1)
    D = scaling_factor * tr(B_train * CCXXCC * B_test' * Σinv)
end


"""
    cvMANOVA(m,eeg)
# Arguments
- `m`: A fitted UnfoldModel with a `LinearModelFitCV` modelfit.
- `eeg`: EEG data used to estimate the noise covariance matrix from the training data (chan x time x trials). Subset if you only want to use e.g. the baseline.
"""
function cvMANOVA(m, eeg; kwargs...)
    @assert isa(modelfit(m), LinearModelFitCV) "cvMANOVA is only implemented for LinearModelFitCV, but got $(typeof(modelfit(m)))"
    cvMANOVA.(Ref(m), Ref(eeg), 1:length(modelfit(m).folds[1]); kwargs...)
end


function cvMANOVA(m, eeg, idx_cv_fold; kwargs...)
    @assert isa(modelfit(m), LinearModelFitCV) "cvMANOVA is only implemented for LinearModelFitCV, but got $(typeof(modelfit(m)))"
    @assert length(modelfit(m).folds) == 1 "cvMANOVA currently does not support multi-event models"

    train_idx = modelfit(m).folds[1][idx_cv_fold].train
    test_idx = modelfit(m).folds[1][idx_cv_fold].test
    X_test = modelmatrix(m)[1][test_idx, :] # needed for cvMANOVA algo


    X_train = modelmatrix(m)[1][train_idx, :] # for noise cov
    Y_train = eeg[:, :, train_idx] # for noise cov


    β_train = modelfit(m).estimate[:, :, :, idx_cv_fold]
    β_test = modelfit(m).info[1].test_estimates[:, :, :, idx_cv_fold]
    cvMANOVA(β_train, β_test, X_train, X_test, Y_train; kwargs...)
end

"""
    cvMANOVA(β_train,β_test,X_train,X_test,Y_train;C, C_test=C,λ=0.05)
A function to perform cross-validated MANOVA (cvMANOVA) on the training and test estimates from a cross validated modelfit
This function calculates the cvMANOVA D statistic for each time point based on the training and test estimates, the design matrices.
It estimates the noise covariance matrix estimated from the given training data.

# Arguments
- `β_train`: The parameter estimates from the training data (chan, time, coeff).
- `β_test`: The parameter estimates from the test data (chan, time, coeff).
- `X_train`: The design matrix for the training data.       
- `X_test`: The design matrix for the test data.
- `Y_train`: The observed data for the training set, used to estimate the noise covariance - subset if you only want to use e.g. the baseline.
# Keyword Arguments
- `C`: The contrast matrix for the training data (no default).
- `C_test`: The contrast matrix for the test data (default: same as C).
- `λ`: Regularization parameter for the noise covariance estimation (default: 0.0 - recommended to be quite small in 2014 paper)
- `temporal_generalization`: calculate all training <-> testing time points, resulting in a matrix of output instead of a vector (the diagonal of the matrix would be the vector). Default off
"""
function cvMANOVA(
    β_train::AbstractArray{T,3},
    β_test::AbstractArray{T,3},
    X_train::AbstractMatrix,
    X_test::AbstractMatrix,
    Y_train::AbstractArray{T,3};
    C::AbstractVector,
    C_test::AbstractVector = C,
    λ = 0.0,
    temporal_generalization = false
) where {T}
    # calculate noise regularized residual error covariance matrix
    Y_train_avg = dropdims(mean(Y_train, dims = 2), dims = 2)
    B_train_avg = X_train \ Y_train_avg'

    @debug size(Y_train_avg) size(X_train) size(B_train_avg)
    Ξ = Y_train_avg' - (X_train * B_train_avg) # "noise" residuals

    Σinv = calculate_Σinv(Ξ; λ)

    # precalculate contrasts & contrast projection
    CC_train = C * pinv(C)
    CC_test = C * pinv(C_test) # no typo, really train * test ^-1

    # unclear why I need this
    csfct = sum(C_test .!= 0) / sum(C .!= 0)
    CCXXCC = CC_train * X_test' * X_test * CC_test / csfct

    # calculate scaling factor
    n_train = size(X_train, 1)
    #n_test = size(X_test, 1)
    n_test = sum(sum(X_test * CC_test .!= 0; dims = 2) .> 0)
    @debug "n_test" size(X_test, 1), n_test


    p = size(β_train, 1) # "voxels", channels here
    fE = n_train - rank(X_test) # residual degree of freedom
    scaling_factor = (fE - p - 1) / n_test


    # calculate D for each time-point
    
    time_idx_train = 1:size(β_test, 2)
    time_idx_test = time_idx_train # assumes length(\beta_test) == length(\beta_train), but that should be given
    if temporal_generalization
        time_idx_train = repeat(time_idx_train, inner = length(time_idx_train)) #[1,2,3,1,2,3,1,2,3]
        time_idx_test = repeat(time_idx_train, outer = length(time_idx_test)) #[1,1,1,2,2,2,3,3,3]
    end

    D = map(
        (t1, t2) -> cvmanova_D(
            CCXXCC,
            Σinv,
            β_train[:, t1, :],
            β_test[:, t2, :];
            scaling_factor,
        ),
        time_idx_train, 
        time_idx_test,
    )
    return temporal_generalization ?  reshape(D, size(β_test,2),size(β_test,2)) : D # return the Vector directly
    
end
