using Revise
using BondGraphs
using LightGraphs
using Test

@testset "BondGraphs.jl" begin
    include("graphfunctions_tests.jl")
    include("construction_tests.jl")
end