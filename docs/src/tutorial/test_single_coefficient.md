````@example eeg
using Random
using CairoMakie
using UnfoldSim
using Unfold
using UnfoldMakie
using Statistics
using ClusterDepth
````

# # Unfold single parameter testing

This is an adaptation of the [ClusterDepth.jl tutorial](https://www.s-ccs.de/ClusterDepth.jl/dev/tutorials/eeg/).

Let's create data from 20 subjects
````@example eeg
data,events = UnfoldSim.predef_eeg(20;return_epoched=true)
times = range(0,step=1/100,length=size(data,2))
````
Fit an UnfoldModel to each subject
````@example eeg
formula = @formula(0 ~ 1 + condition)
models = map((d, ev) -> (fit(UnfoldModel, formula, DataFrame(ev), d, ), ev.subject[1]),
    eachslice(data; dims=3),
    groupby(events, :subject))
````

now we can inspect the data easily, and extract the face-effect

````@example eeg
function add_subject!(df, s)
    df[!, :subject] .= s
    return df
end
allEffects = map((x) -> (effects(Dict(:condition => ["car", "face"]), x[1]), x[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)

plot_erp(allEffects; mapping=(color=:condition, group=:subject))
````

extract the face-coefficient from the linear model

````@example eeg
allCoefs = map(m -> (coeftable(m[1]), m[2]) |> (x) -> add_subject!(x[1], x[2]), models) |> e -> reduce(vcat, e)
plot_erp(allCoefs; mapping=(group=:subject, col=:coefname))
````

let's unstack the tidy-coef table into a matrix and put it to clusterdepth for clusterpermutation testing

````@example eeg
faceCoefs = allCoefs |> x -> subset(x, :coefname => x -> x .== "condition: face")
erpMatrix = unstack(faceCoefs, :subject, :time, :estimate) |> x -> Matrix(x[:, 2:end])' |> collect
summary(erpMatrix)
````

## Clusterdepth

````@example eeg
pvals = clusterdepth(erpMatrix; τ=quantile(TDist(n_subjects - 1), 0.95), nperm=5000);
nothing #hide
````

well - that was fast, less than a second for a cluster permutation test. not bad at all!

## Plotting
Some plotting, and we add the identified cluster

first calculate the ERP

````@example eeg
faceERP = groupby(faceCoefs, [:time, :coefname]) |>
          x -> combine(x, :estimate => mean => :estimate,
    :estimate => std => :stderror);
nothing #hide
````

put the pvalues into a nicer format

````@example eeg
pvalDF = ClusterDepth.cluster(pvals .<= 0.05) |> x -> DataFrame(:from => x[1] ./ 250, :to => (x[1] .+ x[2]) ./ 250, :coefname => "condition: face")
plot_erp(faceERP; stderror=true, pvalue=pvalDF)
````

Looks good to me! We identified the cluster :-)