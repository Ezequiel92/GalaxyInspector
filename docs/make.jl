push!(LOAD_PATH, "./src/")
using Documenter, GalaxyInspector
using DocumenterTools: Themes

# True if it is being deployed, false if it is being compile locally
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# Compile themes
Themes.compile(
    joinpath(@__DIR__, "src/assets/light.scss"),
    joinpath(@__DIR__, "src/assets/themes/documenter-light.css"),
)
Themes.compile(
    joinpath(@__DIR__, "src/assets/dark.scss"),
    joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"),
)

# Compile documentation
makedocs(
    sitename="GalaxyInspector.jl",
    authors="Ezequiel Lozano",
    format   = Documenter.HTML(
        prettyurls=CI,
        assets=["./assets/GI-docs.css"],
        warn_outdated=true,
        collapselevel=1,
        edit_link="main",
    ),
    modules=[GalaxyInspector],
    pages=[
        "Introduction" => "intro.md",
        "API"          => Any[
            "api/convenience.md",
            "api/pipelines.md",
            "api/data_acquisition.md",
            "api/data_analysis.md",
            "api/post_processing.md",
            "api/arepo_utilities.md",
            "api/general_utilities.md",
            "api/constants.md",
        ],
        "Index"        => "index.md",
    ],
)

# Deploy the documentation to GitHub pages
if CI
    deploydocs(
        repo="github.com/Ezequiel92/GalaxyInspector.git",
        devbranch="main",
        versions=["dev" => "dev"],
    )
end
