using Test
using BondGraphs
using BondGraphs: is_connected, components
using BondGraphs: resistor, capacitor, inductor, voltagesource, currentsource, transformer, KCL, KVL
using BondGraphs: chemicalspecies, reaction, chemostat

using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using Graphs
using MetaGraphsNext
using DataStructures
using Catalyst
using DifferentialEquations: Rosenbrock23

function RCI()
    @named c = capacitor()
    @named r = resistor()
    @named i = inductor()
    @named v = voltagesource()
    @named kvl = KVL()
    b1 = Bond(v, kvl)
    b2 = Bond(kvl, c)
    b3 = Bond(kvl, r)
    b4 = Bond(kvl, i)
    BondGraph([c, r, i, v, kvl], [b1, b2, b3, b4], name=:RCI)
end

@testset begin
    @testset "Construction Tests" include("./construction_tests.jl")
    @testset "Equation Tests" include("./equation_tests.jl")
    @testset "Catalyst Tests" include("./catalyst_tests.jl")
    # @testset "Simulation Tests" include("./simulation_tests.jl")
    # @testset "Module Tests" include("./module_tests.jl")
    @testset "Causality Tests" include("./causality_tests.jl")
end


#------------------
# scratchpad
# TODO continue through tests

using BondGraphs: resistor, capacitor, inductor, KVL
using BondGraphs: voltagesource, currentsource, transformer, KCL

@named cap = capacitor()

voltagesource(; E=t->1, name=:v)

transformer(2)

@named kvl = KVL()
