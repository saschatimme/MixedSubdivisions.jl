using Documenter, MixedSubdivisions

makedocs(
    sitename = "MixedSubdivisions.jl",
    pages = [
        "MixedSubdivisions" => "index.md",
    ]
)

deploydocs(
    repo   = "github.com/JuliaHomotopyContinuation/MixedSubdivisions.jl.git",
)
