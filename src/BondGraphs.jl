module BondGraphs

import LightGraphs; const lg = LightGraphs
import Base: RefValue, eltype, show, in

export AbstractNode, Component, Junction, Bond, BondGraph,
vertex, set_vertex!,
add_node!, remove_node!, connect!, disconnect!, swap!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")

end
