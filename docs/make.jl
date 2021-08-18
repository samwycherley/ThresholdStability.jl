using ThresholdStability
using Documenter

DocMeta.setdocmeta!(ThresholdStability, :DocTestSetup, :(using ThresholdStability); recursive=true)

makedocs(;
    modules=[ThresholdStability],
    authors="Sam Wycherley",
    repo="https://github.com/samwycherley/ThresholdStability.jl/blob/{commit}{path}#{line}",
    sitename="ThresholdStability.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://samwycherley.github.io/ThresholdStability.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Stability" => "stability.md"
        "CKSVAR" => "cksvar.md"
    ],
)

deploydocs(;
    repo="github.com/samwycherley/ThresholdStability.jl",
)
