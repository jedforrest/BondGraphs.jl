# Introduction

## Overview
[BondGraphs.jl](https://github.com/jedforrest/BondGraphs.jl) is a Julia implementation of the bond graph framework. Bond graphs are an energy-based modelling framework that describe energy flow through a physical system, and are especially useful for modelling multi-scale or multi-physical systems.[^1]

[^1]: Gawthrop and Bevan, _Bond-graph modeling_ (2007)

This package constructs a symbolic graph model of a physical system, which can then be converted into a system of differential equations. BondGraphs.jl includes specific methods that interact with the wider Julia modelling ecosystem, including ModelingToolkit.jl, Catalyst.jl, and Plots.jl.

```@contents
Pages   = [
    "gettingstarted.md",
    "examples.md",
    "api.md"
]
Depth = 1
```

## Installation
```julia
using Pkg; Pkg.add("BondGraphs")
using BondGraphs
```

## Tutorials
For tutorials and examples, refer to the [Examples](@ref) page. For interactive Jupyter Notebook tutorials, refer to [BondGraphsTutorials](https://github.com/jedforrest/BondGraphsTutorials) on GitHub.
