using Revise
using BondGraphs
using LightGraphs
using Catalyst
using Test

@testset "BondGraphs.jl" begin
    include("graphfunctions_tests.jl")
    include("construction_tests.jl")
    include("conversion_tests.jl")
end