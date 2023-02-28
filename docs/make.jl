using Documenter, BondGraphs

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
