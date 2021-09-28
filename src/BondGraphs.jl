module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in, ==

using StaticArrays
using Catalyst

export AbstractNode, Component, Junction, Port, Bond, BondGraph,
vertex, set_vertex!, freeports, numports, srcnode, dstnode, getnodes, getbonds,
add_node!, remove_node!, connect!, disconnect!, swap!, insert_node!, merge!

include("basetypes.jl")
include("basetypefunctions.jl")
include("graphfunctions.jl")
include("construction.jl")
include("conversion.jl")

end
