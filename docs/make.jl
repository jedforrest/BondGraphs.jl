using Documenter, BondGraphs, Graphs, Plots, ModelingToolkit, Catalyst

makedocs(
    sitename = "BondGraphs.jl",
    modules = [BondGraphs],
    pages = [
        "index.md",
        # "background.md",
        "gettingstarted.md",
        "examples.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/jedforrest/BondGraphs.jl.git",
)
