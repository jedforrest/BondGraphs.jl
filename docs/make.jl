using Documenter, BondGraphs

makedocs(
    sitename = "BondGraphs.jl",
    modules = [BondGraphs],
    pages = [
        "index.md",
        # "background.md",
        "gettingstarted.md",
        "userguide.md",
        "examples.md",
        "API Reference" => "api.md"
    ]
)
