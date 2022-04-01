push!(LOAD_PATH, "../src/")
using Documenter, GadgetInspector
using DocumenterTools: Themes

# True if being deployed, false if being compile locally
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# Compile themes
Themes.compile(
    joinpath(@__DIR__, "src/assets/geostats-light.scss"), 
    joinpath(@__DIR__, "src/assets/themes/documenter-light.css"),
)
Themes.compile(
    joinpath(@__DIR__, "src/assets/geostats-dark.scss"), 
    joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"),
)

# Compile documentation
makedocs(
    sitename = "GadgetInspector.jl",
    authors = "Ezequiel Lozano",
    format = Documenter.HTML(
        assets=[
            asset( 
				"https://fonts.googleapis.com/css?family=Inconsolata&display=swap",
                class = :css,
            )
        ],
        prettyurls = CI,
        edit_link = "main",
    ),
    modules = [GadgetInspector],
    pages = [
        "Introduction"                  => "intro.md",
        "Pipelines"                     => "pipelines.md",
		"Data acquisition functions"    => "data_acquisition.md",
		"Data analysis functions"       => "data_analysis.md",
		"Post processing functions"     => "post_processing.md",
		"Utilities"                     => "utilities.md",
		"Constants"                     => "constants.md",
        "Index"                         => "index.md",
    ],
)

# Deploy documentation to GitHub pages
if CI
    deploydocs(
        repo =  "github.com/Ezequiel92/GadgetInspector.git",
        devbranch = "main",
        versions = ["dev" => "dev"],
    ) 
end
