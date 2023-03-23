using Test
using BondGraphs
using Graphs
using ModelingToolkit
using DifferentialEquations: Rosenbrock23
using Catalyst
using RecipesBase

@testset "BondGraphs.jl" begin
    @testset "Graph functions" begin include("graphfunctions_tests.jl") end
    @testset "Construction" begin include("construction_tests.jl") end
    @testset "Equations" begin include("equation_tests.jl") end
    @testset "Simulations" begin include("simulation_tests.jl") end
    @testset "Catalyst" begin include("catalyst_tests.jl") end
    @testset "Miscellaneous" begin include("misc_tests.jl") end
    @testset "Modules" begin include("module_tests.jl") end
end
