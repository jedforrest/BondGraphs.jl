module BondGraphs

import LightGraphs; const lg = LightGraphs
import Base: RefValue, eltype, show

export AbstractNode, Component, Junction, Bond, BondGraph,
vertex, set_vertex!,
find_index, add_nodes!, connect!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")

end
