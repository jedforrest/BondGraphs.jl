module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in, ==

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils
using SymbolicUtils.Rewriters
using DataStructures
using Setfield

export AbstractNode, Component, Junction, Port, Bond, BondGraph,
EqualEffort, EqualFlow,
vertex, set_vertex!, freeports, numports, srcnode, dstnode, 
cr, params, state_vars, set_param!, set_initial_value!, 
default_value, equations, simulate,
new, add_node!, remove_node!, connect!, disconnect!, swap!

include("basetypes.jl")
include("graphfunctions.jl")
include("construction.jl")
include("components.jl")
include("equations.jl")

end
