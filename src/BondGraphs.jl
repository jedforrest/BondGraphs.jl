module BondGraphs

import Graphs as g
import Base: RefValue, eltype, show, in, iterate, ==, getproperty
import ModelingToolkit: parameters, states, equations # Importing means names can be reused, but may be confusing

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst

export AbstractNode, Component, Junction, EqualEffort, EqualFlow,
Port, Bond, BondGraph, BondGraphNode,

type, name, freeports, numports, weights, vertex, set_vertex!,
parameters, states, defaults, constitutive_relations,
get_parameter, set_parameter!, get_initial_value, set_initial_value!,

srcnode, dstnode, nodes, bonds, getnodes, getbonds,

add_node!, remove_node!, connect!, disconnect!, 
swap!, insert_node!, merge_nodes!, simplify_junctions!,

simulate, addlibrary!

# Component libraries
include("libraries/biochemical.jl")
include("libraries/standard.jl")
include("libraries/libraryfunctions.jl")

include("basetypes/AbstractNode.jl")
include("basetypes/Bond.jl")
include("basetypes/BondGraph.jl")

include("graphfunctions.jl")
include("construction.jl")
include("equations.jl")
include("conversion.jl")

end
