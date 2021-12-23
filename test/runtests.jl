using Revise
using BondGraphs
using BondGraphs.BiochemicalLibrary
using Graphs
using ModelingToolkit
using Catalyst
using Test

@testset "BondGraphs.jl" begin
    @testset "Graph functions" begin include("graphfunctions_tests.jl") end
    @testset "Construction" begin include("construction_tests.jl") end
    @testset "Equations" begin include("equation_tests.jl") end
    # @testset "Simulations" begin include("simulation_tests.jl") end
    @testset "Conversions" begin include("conversion_tests.jl") end
end