module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in

using StaticArrays
using ModelingToolkit
using SymbolicUtils

export AbstractNode, Component, Junction, Port, Bond, BondGraph,
vertex, set_vertex!, freeports, numports, srcnode, dstnode,
new, add_node!, remove_node!, connect!, disconnect!, swap!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")
include("components.jl")

end
