using Test
using BondGraphs
using BondGraphs: is_connected
using BondGraphs: resistor, capacitor, inductor, voltagesource, KVL
using BondGraphs: chemicalspecies, reaction, chemostat

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Graphs
using MetaGraphsNext
using DataStructures
# using DifferentialEquations: Rosenbrock23
# using Catalyst
# using RecipesBase

function RCI()
    @named r = resistor()
    @named c = capacitor()
    @named i = inductor()
    @named v = voltagesource()
    @named kvl = KVL()
    b2 = Bond(r, kvl)
    b1 = Bond(c, kvl)
    b3 = Bond(kvl, i)
    b4 = Bond(kvl, v)
    BondGraph([c, r, i, v, kvl], [b1, b2, b3, b4], name=:RCI)
end

@testset begin
    @testset "Construction Tests" include("./construction_tests.jl")
    @testset "Equation Tests" include("./equation_tests.jl")
    # @testset "Catalyst Tests" include("./catalyst_tests.jl")
    # @testset "Simulation Tests" include("./simulation_tests.jl")
    # @testset "Module Tests" include("./module_tests.jl")
end
