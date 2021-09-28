using Revise
using BondGraphs
using LightGraphs
using ModelingToolkit
using Test

@testset "BondGraphs.jl" begin
    include("graphfunctions_tests.jl")
    include("construction_tests.jl")
    include("equation_tests.jl")
end