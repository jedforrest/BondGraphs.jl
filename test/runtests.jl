using Test
using TestItems
using TestItemRunner
using BondGraphs
using Graphs
using ModelingToolkit
using DifferentialEquations: Rosenbrock23
using Catalyst
using RecipesBase

@run_package_tests

@testset "Graph functions" begin include("graphfunctions_tests.jl") end
# @testitem "Construction" begin include("construction_tests.jl") end
# @testitem "Equations" begin include("equation_tests.jl") end
# @testitem "Simulations" begin include("simulation_tests.jl") end
# @testitem "Catalyst" begin include("catalyst_tests.jl") end
# @testitem "Miscellaneous" begin include("misc_tests.jl") end
# @testitem "Modules" begin include("module_tests.jl") end
