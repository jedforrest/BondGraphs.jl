module BondGraphs

import Graphs as g
import Base: RefValue, eltype, show, in, iterate, ==, getproperty
# Importing means names can be reused, but may be confusing later
import ModelingToolkit: parameters, states, equations, controls

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst
using RecipesBase, GraphRecipes

export AbstractNode, Component, Junction, EqualEffort, EqualFlow,
SourceSensor, Port, Bond, BondGraph, BondGraphNode,

type, name, freeports, numports, weights, vertex, set_vertex!,
parameters, states, controls, constitutive_relations,
get_default, set_default!,

srcnode, dstnode, nodes, bonds, components, junctions, getnodes, getbonds,

add_node!, remove_node!, connect!, disconnect!, 
swap!, insert_node!, merge_nodes!, simplify_junctions!, expose,

simulate, addlibrary!, description

# Component libraries
include("libraries/biochemical.jl")
include("libraries/standard.jl")
include("libraries/libraryfunctions.jl")

# Types used by BondGraphs
include("basetypes/AbstractNode.jl")
include("basetypes/Bond.jl")
include("basetypes/BondGraph.jl")

# Core functionality
include("graphfunctions.jl")
include("construction.jl")
include("equations.jl")
include("conversion.jl")
include("plotrecipes.jl")

end
