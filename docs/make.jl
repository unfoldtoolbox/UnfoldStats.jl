using UnfoldStats
using Documenter

DocMeta.setdocmeta!(UnfoldStats, :DocTestSetup, :(using UnfoldStats); recursive=true)

makedocs(;
    modules=[UnfoldStats],
    authors="Judith Schepers, Benedikt V. Ehinger",
    repo="https://github.com/unfoldtoolbox/UnfoldStats.jl/blob/{commit}{path}#{line}",
    sitename="UnfoldStats.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://unfoldtoolbox.github.io/UnfoldStats.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/unfoldtoolbox/UnfoldStats.jl",
    devbranch="main",
)
