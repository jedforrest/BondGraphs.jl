# BondGraphs.jl

![CI](https://github.com/jedforrest/BondGraphs.jl/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/jedforrest/BondGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jedforrest/BondGraphs.jl)

BondGraphs.jl is a Julia implementation of the bond graph framework, built upon existing packages [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl), and [Catalyst.jl](https://github.com/SciML/Catalyst.jl). This package constructs a symbolic graph model of a physical system, which can then be converted into a system of differential equations.

For tutorials and code demonstrations in Jupyter Notebooks, see [BondGraphsTutorials](https://github.com/jedforrest/BondGraphsTutorials).

## Installation
```julia
using Pkg; Pkg.add(url="https://github.com/jedforrest/BondGraphs.jl")
using BondGraphs
```

## Bond graphs
Bond graphs are an energy-based modelling framework that describe the rate of energy flow moving through system components. By construction, bond graph models enforce physical and thermodynamic constraints, guaranteeing compatibility with other physical models. This framework has been applied to mechanical, electrical, chemical, and biological systems, and is even capable of modelling complex multi-physics systems.




