module BondGraphs

import LightGraphs as lg
import Base: RefValue, eltype, show, in, ==, getproperty

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils
using SymbolicUtils.Rewriters
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

include("basetypes.jl")
#include("basetypefunctions.jl")
include("graphfunctions.jl")
include("construction.jl")
include("components.jl")
include("equations.jl")
include("conversion.jl")

end
