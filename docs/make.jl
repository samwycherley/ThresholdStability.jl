using ThresholdStability
using Documenter

makedocs(root="/",
    source="src",
    build="build",
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
    ],
)

deploydocs(;
    repo="github.com/samwycherley/ThresholdStability.jl",
)
