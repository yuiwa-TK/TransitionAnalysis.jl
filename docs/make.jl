using TransitionAnalysis
using Documenter

DocMeta.setdocmeta!(TransitionAnalysis, :DocTestSetup, :(using TransitionAnalysis); recursive=true)

makedocs(;
    modules=[TransitionAnalysis],
    authors="Yuta Iwatani",
    repo="https://github.com/iwatani/TransitionAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="TransitionAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iwatani.github.io/TransitionAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iwatani/TransitionAnalysis.jl",
    devbranch="main",
)
