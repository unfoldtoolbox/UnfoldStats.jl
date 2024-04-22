using UnfoldStats
using Documenter
using Literate
using Glob

GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__, "literate")

# Check whether the "generated" directory already exists, otherwise create it
if !isdir(GENERATED)
    mkdir(GENERATED)
    @info "`generated` directory has been created."
end

for subfolder ∈ ["explanation", "how-to", "tutorial", "reference"]
    local SOURCE_FILES = Glob.glob(subfolder * "/*.jl", SOURCE)
    foreach(fn -> Literate.markdown(fn, GENERATED * "/" * subfolder), SOURCE_FILES)

end

DocMeta.setdocmeta!(UnfoldStats, :DocTestSetup, :(using UnfoldStats); recursive = true)

makedocs(;
    modules = [UnfoldStats],
    authors = "Judith Schepers, Benedikt V. Ehinger",
    repo = "https://github.com/unfoldtoolbox/UnfoldStats.jl/blob/{commit}{path}#{line}",
    sitename = "UnfoldStats.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://unfoldtoolbox.github.io/UnfoldStats.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" =>
            ["Two-stage EEG analysis" => "generated/tutorial/two_stage_analysis.md"],
        ["LMM EEG + ClusterDepth" => "generated/tutorial/lmm_clusterdepth.md"],
    ],
)

deploydocs(; repo = "github.com/unfoldtoolbox/UnfoldStats.jl", devbranch = "main")
