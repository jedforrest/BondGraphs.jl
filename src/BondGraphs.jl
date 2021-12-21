module BondGraphs

import Graphs as g
import Base: RefValue, eltype, show, in, iterate, ==, getproperty

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst

export AbstractNode, Component, Junction, EqualEffort, EqualFlow,
Port, Bond, BondGraph, BondGraphNode,

type, name, portconnections, portweights, numports, vertex, set_vertex!,
params, state_vars, equations,

srcnode, dstnode, nodes, bonds, getnodes, getbonds,

add_node!, remove_node!, connect!, disconnect!, 
swap!, insert_node!, merge_nodes!, simplify_junctions!

# cr, set_param!, set_initial_value!, simulate

# Component libraries
include("libraries/biochemical.jl")
include("libraries/standard.jl")
using .StandardLibrary

include("basetypes/AbstractNode.jl")
include("basetypes/Bond.jl")
include("basetypes/BondGraph.jl")

include("graphfunctions.jl")
include("construction.jl")
# include("equations.jl")
include("conversion.jl")

end
