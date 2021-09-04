module BondGraphs

import LightGraphs; const lg = LightGraphs
import Base: RefValue, eltype, show, in

using StaticArrays

export AbstractNode, Component, Junction, Port, Bond, BondGraph,
vertex, set_vertex!, freeports, numports, srcnode, dstnode,
add_node!, remove_node!, connect!, disconnect!, swap!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")

end
