__precompile__(false)

module BondGraphs

import Base: RefValue, eltype, show, in, iterate, ==, getproperty, setproperty!, getindex, setindex!
import Graphs:
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
# import ModelingToolkitBase: parameters, equations, controls

using Graphs: SimpleGraph, SimpleDiGraph, SimpleEdge
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using DifferentialEquations
using SymbolicUtils, SymbolicUtils.Rewriters
using OrderedCollections
using Setfield
using Catalyst
using RecipesBase, GraphRecipes
using Latexify
using Graphs
using GraphMakie, GraphMakie.NetworkLayout
using MetaGraphsNext
using DataStructures


const Effort = ModelingToolkitBase.Equality

# Structures
include("structures/AbstractElement.jl")
include("structures/Port.jl")
include("structures/Component.jl")
include("structures/Bond.jl")
include("structures/BondGraph.jl")
include("structures/BondGraphElement.jl")

# Core functionality
# include("graphfunctions.jl")
include("construction.jl")
# include("systems.jl")
include("catalyst.jl")
include("plotrecipes.jl")

# Component libraries
include("libraries/ComponentLibrary.jl")


export Component,
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
    srccomp,
    dstcomp,
    srclabel,
    dstlabel,
    comps,
    bonds,
    components,
    junctions,
    getcomps,
    getbonds,
    add_comp!,
    remove_comp!,
    connect!,
    disconnect!,
    swap!,
    insert_comp!,
    merge_comps!,
    simplify_junctions!,
    expose,
    simulate,
    addlibrary!,
    description,
    bgplot,
    # NEW
    Port,
    Effort,
    Compliance,
    Inertance,
    Resistance,
    EffortSource,
    FlowSource,
    SourceSensor,
    Transformer,
    Gyrator,
    EqualEffort,
    EqualFlow,
    ZeroJunction,
    OneJunction,
    power_variables,
    elementtype,
    efforts,
    flows,
    states,
    system,
    is_connected,
    components,
    connect!,
    effort,
    flow,
    elements,
    junctions,
    componentnames,
    inneighbor_comps,
    outneighbor_comps,
    compile_system!,
    sub_defaults,
    BondGraphElement,
    assign_causality

end
