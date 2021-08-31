using ThresholdStability
using Documenter

# DocMeta.setdocmeta!(ThresholdStability, :DocTestSetup, :(using ThresholdStability); recursive=true)

makedocs(
    authors="Sam Wycherley",
    sitename="ThresholdStability.jl",
    format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true"
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
