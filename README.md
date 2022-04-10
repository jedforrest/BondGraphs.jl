# BondGraphs

[![Tests status](https://ci.appveyor.com/api/projects/status/github/jedforrest/bondgraphs?svg=true)](https://ci.appveyor.com/project/jedforrest/bondgraphs)
[![codecov](https://codecov.io/gh/jedforrest/BondGraphs/branch/master/graph/badge.svg)](https://codecov.io/gh/jedforrest/BondGraphs)

Bond graphs are an energy-based modelling framework that describe the rate of energy flow moving through system components. By construction, bond graph models enforce physical and thermodynamic constraints, guaranteeing compatibility with other physical models. This framework has been applied to mechanical, electrical, chemical, and biological systems, and is even capable of modelling complex multi-physics systems.

BondGraphs.jl is a pure Julia implementation of the bond graph framework, built upon existing packages [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl). This package constructs a symbolic graph model of a physical system, which can then be converted into a system of differential equations. The models are composable and hierarchical, and compatible with the Julia modelling ecosystem.

For tutorials and code demonstrations in Jupyter Notebooks, see [BondGraphsTutorials](https://github.com/jedforrest/BondGraphsTutorials).
