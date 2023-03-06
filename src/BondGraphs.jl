module BondGraphs

import Graphs as g
import Base: RefValue, eltype, show, in, iterate, ==, getproperty, setproperty!
import ModelingToolkit: parameters, states, equations, controls

using StaticArrays
using ModelingToolkit
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst
using RecipesBase, GraphRecipes
using Latexify

export AbstractNode, Component, Junction, EqualEffort, EqualFlow,
SourceSensor, Port, Bond, BondGraph, BondGraphNode,

type, name, id, freeports, numports, weights, vertex, set_vertex!,
parameters, globals, states, controls, all_variables, constitutive_relations,
has_controls,

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
include("catalyst.jl")
include("plotrecipes.jl")

end
