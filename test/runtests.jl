using Revise
using BondGraphs
using LightGraphs
using Catalyst
using Test

@testset "BondGraphs.jl" begin
    @testset "Basic Graph Functionality" begin
        include("graphfunctions_tests.jl")
    end
    @testset "Bond Graph Construction" begin
        include("construction_tests.jl")
    end
    @testset "Reaction Network Conversion" begin
        include("conversion_tests.jl")
    end
end