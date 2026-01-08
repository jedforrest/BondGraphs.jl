__precompile__(false)

module BondGraphs

import Base: RefValue, eltype, show, in, iterate, ==, getproperty, setproperty!, getindex
import Graphs:
               SimpleGraph,
               SimpleDiGraph,
               edgetype,
               edges,
               ne,
               has_edge,
               vertices,
               nv,
               has_vertex,
               inneighbors,
               outneighbors,
               all_neighbors,
               is_directed,
               zero,
               src,
               dst,
               add_vertex!,
               rem_vertex!,
               add_edge!,
               rem_edge!
import ModelingToolkit: parameters, equations, controls

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst
using RecipesBase, GraphRecipes
using Latexify
using Graphs
using GraphMakie, GraphMakie.NetworkLayout

export AbstractNode,
       Component,
    #    Junction,
       EqualEffort,
       EqualFlow,
       SourceSensor,
       Bond,
       BondGraph,
    #    BondGraphNode,
       type,
       name,
       id,
       ports,
       numports,
       isconnected,
       weights,
       vertex,
       set_vertex!,
       parameters,
       globals,
       states,
       controls,
       all_variables,
       constitutive_relations,
       has_controls,
       srcnode,
       dstnode,
       srclabel,
       dstlabel,
       nodes,
       bonds,
       components,
       junctions,
       getnodes,
       getbonds,
       add_node!,
       remove_node!,
       connect!,
       disconnect!,
       swap!,
       insert_node!,
       merge_nodes!,
       simplify_junctions!,
       expose,
       simulate,
       addlibrary!,
       description,
       bgplot,
       # NEW
       Port,
       Effort,
       StaticStorageElement,
       DynamicStorageElement,
       DissipatorElement,
       EffortSource,
       FlowSource,
       SourceSensor,
       Transformer,
       Gyrator,
       EqualEffort,
       EqualFlow,
       power_variables,
       elementtype,
       efforts,
       flows,
       states,
       system,
       is_connected,
       connect!,
       effort,
       flow,
       elements,
       junctions

# Component libraries
include("libraries/biochemical.jl")
include("libraries/standard.jl")
include("libraries/libraryfunctions.jl")

# Structures
include("structures/Port.jl")
include("structures/AbstractElement.jl")
include("structures/AbstractNode.jl")
include("structures/Bond.jl")
include("structures/BondGraph.jl")

# Core functionality
include("graphfunctions.jl")
# include("construction.jl")
include("systems.jl")
include("catalyst.jl")
include("plotrecipes.jl")

# NEW
include("libraries/electrical.jl")

end
