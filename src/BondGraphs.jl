module BondGraphs

import LightGraphs; const lg = LightGraphs
import Base: eltype, show

export AbstractNode, Component, Junction, Bond, BondGraph

include("basetypes.jl")
include("graphfunctions.jl")

end
