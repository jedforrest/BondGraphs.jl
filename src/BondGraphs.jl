module BondGraphs

import Graphs as g
# import ModelingToolkit as mtk
import Base: RefValue, eltype, show, in, iterate, ==, getproperty

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using DataStructures
using Setfield
using Catalyst

export AbstractNode, Component, Junction, Port, Bond, BondGraph, BondGraphNode,
EqualEffort, EqualFlow,

type, name, vertex, set_vertex!, freeports, numports, bondgraph, nodes, bonds,
srcnode, dstnode, getnodes, getbonds,

new, add_node!, remove_node!, connect!, disconnect!, 
swap!, insert_node!, merge_nodes!, simplify_junctions!,

cr, params, state_vars, set_param!, set_initial_value!, 
default_value, equations, simulate

include("basetypes/AbstractNode.jl")
include("basetypes/Bond.jl")
include("basetypes/BondGraph.jl")
#include("basetypefunctions.jl")
include("graphfunctions.jl")
include("construction.jl")
include("components.jl")
include("equations.jl")
include("conversion.jl")

end
