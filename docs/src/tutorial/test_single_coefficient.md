

# Unfold single parameter cluster permutation testing
This is an adaptation of the [ClusterDepth.jl tutorial](https://www.s-ccs.de/ClusterDepth.jl/dev/tutorials/eeg/).

We simulate a 1x2 design and perform cluster permutation testing via the `ClusterDepth.jl` package

## Simulate data
Let's create data from 20 subjects

 ```@raw html
 <details>
 <summary>Click to expand the simulation details</summary>
 ```
```@example eeg
using Random
using CairoMakie
using UnfoldSim
using Unfold
using UnfoldMakie
using Statistics
using ClusterDepth
using DataFrames
using Distributions
```



```@example eeg
n_subjects = 20
data,events = UnfoldSim.predef_eeg(n_subjects;return_epoched=true)
times = range(0,step=1/100,length=size(data,1))
```

 ```@raw html
 </details >
 ```

## Fitting regression Models
Fit an UnfoldModel to each subject:
```@example eeg
formula = @formula(0 ~ 1 + condition)
models = map((d, ev) -> (fit(UnfoldModel, formula, DataFrame(ev), collect(d), times), ev.subject[1]),
    eachslice(data; dims=3),
    groupby(events, :subject))
nothing #hide

```

now we can inspect the data easily, and extract the face-effect

```@example eeg
function add_subject!(df, s)
    df[!, :subject] .= s
    return df
end
allEffects = map((x) -> (effects(Dict(:condition => ["car", "face"]), x[1]), x[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)

plot_erp(allEffects; mapping=(color=:condition, group=:subject))
```

extract the face-coefficient from the linear model

```@example eeg
allCoefs = map(m -> (coeftable(m[1]), m[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)
plot_erp(allCoefs; mapping=(group=:subject, col=:coefname))
```

let's unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing

```@example eeg
faceCoefs = allCoefs |> x -> subset(x, :coefname => x -> x .== "condition: face")
erpMatrix = unstack(faceCoefs, :subject, :time, :estimate) |> x -> Matrix(x[:, 2:end])' |> collect
summary(erpMatrix)
```

## Clusterdepth

```@example eeg

pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95), nperm=5000);
nothing #hide
```

well - that was fast, less than a second for a cluster permutation test. not bad at all!

## Plotting
Some plotting, and we add the identified cluster

first calculate the ERP

```@example eeg
faceERP = groupby(faceCoefs, [:time, :coefname]) |>
          x -> combine(x, :estimate => mean => :estimate,
    :estimate => std => :stderror);
nothing #hide
```

put the significance into a dataframe-form

```@example eeg
significant_regions_df = ClusterDepth.cluster(pvals .<= 0.05) |> x -> DataFrame(:from => times[x[1]], :to => times[x[1] .+ x[2]] , :coefname => "condition: face")
```

and plot it!
```@example eeg
plot_erp(faceERP; stderror=true, significance=significant_regions_df)
```

Looks good to me! We identified the cluster :-)