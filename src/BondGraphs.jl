module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in, getproperty, ==

using StaticArrays
using Catalyst

export AbstractNode, Component, Junction, Port, Bond, BondGraph, BondGraphNode,
type, name, vertex, set_vertex!, freeports, numports, bondgraph, nodes, bonds,
srcnode, dstnode, getnodes, getbonds,
add_node!, remove_node!, connect!, disconnect!, 
swap!, insert_node!, merge_nodes!, simplify_junctions!

include("basetypes.jl")
include("basetypefunctions.jl")
include("graphfunctions.jl")
include("construction.jl")
include("conversion.jl")

end
