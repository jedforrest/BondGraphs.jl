module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in

using StaticArrays
using ModelingToolkit
using SymbolicUtils
using SymbolicUtils.Rewriters
using DataStructures

export AbstractNode, Component, Junction, Port, Bond, BondGraph,
EqualEffort, EqualFlow,
vertex, set_vertex!, freeports, numports, srcnode, dstnode, 
equations, params, state_vars, control_space, bond_space,
new, add_node!, remove_node!, connect!, disconnect!, swap!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")
include("components.jl")
include("equations.jl")
include("mappings.jl")

end
