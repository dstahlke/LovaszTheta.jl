using Documenter, LovaszTheta

DocMeta.setdocmeta!(LovaszTheta, :DocTestSetup, :( using LovaszTheta; using LinearAlgebra ); recursive=true)

makedocs(
    sitename="LovaszTheta.jl",
    modules=[LovaszTheta],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(
    repo = "github.com/dstahlke/LovaszTheta.jl.git",
    devbranch = "master",
    branch = "gh-pages",
)
