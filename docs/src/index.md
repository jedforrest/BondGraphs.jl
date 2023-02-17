# Introduction

## Overview
[BondGraphs.jl](https://github.com/jedforrest/BondGraphs.jl) is a Julia implementation of the bond graph framework. This package constructs a symbolic graph model of a physical system, which can then be converted into a system of differential equations. Bond graphs are compatible other packages including [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl), and [Catalyst.jl](https://github.com/SciML/Catalyst.jl).

```@contents
Pages   = [
    "gettingstarted.md",
    "userguide.md",
    "examples.md",
    "api.md"
]
Depth = 1
```

## Installation
```julia
using Pkg; Pkg.add(url="https://github.com/jedforrest/BondGraphs.jl")
using BondGraphs
```

## Tutorials
For tutorials and examples, refer to the [Examples](@ref) page. For interactive Jupyter Notebook tutorials, refer to [BondGraphsTutorials](https://github.com/jedforrest/BondGraphsTutorials) on GitHub.
