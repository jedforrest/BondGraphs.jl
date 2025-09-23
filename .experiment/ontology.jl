# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

# Maybe 'component' should be 'element' which includes functionality for
# components and junctions. Then dispatch on the parametric type:
# - StaticStorage
# - DynamicStorage
# - Dissipator
# - Junction
# - Transformer etc.

# Doing it this way means we can get the benefit of dispatch
# while easily able to extend to new sub types for custom components
# e.g. can define a chemical energy store:
#   Ce <: StaticStorage

# These abstract types could define specific components stored in the standard library
# e.g. electric capacitor <: static storage
# these define a fixed definition that is reusable across a bond graph definition
# then the Component struct is specifically for repeated *instances* within a single bond graph model
# TODO See https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects

# TODO I should think of the model-building side of this package as generating MTK models
# from a graph description or interface - I don't need to reinvent the wheel when it comes
# to definining ports and components
import Base: show, size
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

############################################################

abstract type BondGraphVertex end

abstract type BondElement <: BondGraphVertex end
abstract type StorageElement <: BondElement end
abstract type SourceElement <: BondElement end

@variables R
cr = (e, f) -> R * f - e

"""`R` component"""
struct DissipatorElement <: BondElement
    cr::Vector{Any}
end
# function (d::DissipatorElement)(e, f)
#     d.cr(e, f)
# end
# @variables e f
# DissipatorElement(e, f)
# returns R*f - e

# TODO
"""`C` component"""
struct StaticStorageElement <: StorageElement
    cr::Vector{Equation}
end
"""`I` component"""
struct DynamicStorageElement <: StorageElement
    cr::Vector{Equation}
end

"""`Se` component"""
struct EffortSource <: SourceElement end
"""`Sf` component"""
struct FlowSource <: SourceElement end
"""`SS` component"""
struct SourceSensor <: SourceElement end

abstract type JunctionStructure <: BondGraphVertex end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

"""`TF` component"""
struct Transformer <: ParametricJunction end
"""`GY` component"""
struct Gyrator <: ParametricJunction end

"""`0`-junction"""
struct EqualEffort <: NonParametricJunction end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction end

# Julia Type trees
# using GraphRecipes, Plots
# default(size=(1000, 1000))
# plot(BondGraphVertexClass, method=:tree, fontsize=10, nodeshape=:ellipse)

# TODO the above structs include the CR

numports(c::BondElement) = length(c.cr)

# CRs map (e,f) -> ϕ
# TODO define for other vertex types
constitutive_relations(c::BondElement) = c.cr

# Used when displaying in a graph.
# TODO these can just be included in the struct definitions above (kwdef)
glyph(::DissipatorElement) = :R
glyph(::StaticStorageElement) = :C
glyph(::DynamicStorageElement) = :I
glyph(::EffortSource) = :Se
glyph(::FlowSource) = :Sf
glyph(::Transformer) = :TF
glyph(::Gyrator) = :GY
glyph(::JunctionStructure) = :J
glyph(::EqualEffort) = :𝟎
glyph(::EqualFlow) = :𝟏

############################################################
mutable struct Port
    parent::Any  # should be Component
    index::Int
    connected::Bool  # Ref value?
    effort::Num
    flow::Num
    function Port(parent::Any, index::Int=1)
        effort = Symbolics.variable(:e, index)
        flow = Symbolics.variable(:f, index)
        new(parent, index, false, effort, flow)
    end
end
is_connected(p::Port) = p.connected
parent(p::Port) = p.parent
effort(p::Port) = p.effort
flow(p::Port) = p.flow

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "p$(port.index) $connection_state")
end

############################################################

# Components now define BG elements and junction structures
struct Component{V<:BondGraphVertex}
    type::V
    name::Symbol
    ports::Vector{Port}
end
function Component(type::BondElement, name)
    newcomp = Component(type, Symbol(name), Port[])
    for i in 1:numports(type)
        push!(newcomp.ports, Port(newcomp, i))
    end
    newcomp
end
function Component(type::JunctionStructure, name)
    Component(type, Symbol(name), Port[])
end

show(io::IO, comp::Component{<:BondElement}) = print(io, "$(glyph(comp.type))::$(comp.name)")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.type))")

hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, comp.ports))
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    port = Port(junc, index)
    push!(junc.ports, port)
    port
end

constitutive_relations(c::Component) = constitutive_relations(c.type)

############################################################

struct Bond
    src::Port
    dst::Port
    function Bond(src::Port, dst::Port)
        src.connected && error("$src already connected")
        dst.connected && error("$dst already connected")
        src.connected = true
        dst.connected = true
        new(src, dst)
    end
end
function Bond(src::Component, dst::Component)
    hasfreeport(src) || error("$src has no free ports")
    hasfreeport(dst) || error("$dst has no free ports")
    srcport = nextfreeport(src)
    dstport = nextfreeport(dst)
    Bond(srcport, dstport)
end

vertices(b::Bond) = parent(b.src), parent(b.dst)

function show(io::IO, b::Bond)
    src, dst = vertices(b)
    print(io, "$src ⇀ $dst")
end

############################################################

struct BondGraph
    name::Symbol
    elements::Vector{Component}
    junctions::Vector{Component}
    bonds::Vector{Bond}
    function BondGraph(name, elements=[], junctions=[], bonds=[])
        new(Symbol(name), elements, junctions, bonds)
    end
end
# 3 argument constructor
function BondGraph(name, elements_junctions, bonds)
    elements = filter(x -> x.type isa BondElement, elements_junctions)
    junctions = filter(x -> x.type isa JunctionStructure, elements_junctions)
    BondGraph(Symbol(name), elements, junctions, bonds)
end
# 2 argument constructor
function BondGraph(name, bonds::Vector{Bond})
    src_dsts = reduce(vcat, collect.(vertices.(bonds)))
    BondGraph(Symbol(name), unique(src_dsts), bonds)
end
components(bg::BondGraph) = [bg.elements; bg.junctions]

# maybe not as the default "show" option (summary print instead)
function show(io::IO, bg::BondGraph)
    print_str = "BondGraph \"$(bg.name)\""
    print_str *= isempty(bg.bonds) ? "" : "\n$(join(bg.bonds,"\n"))"
    print(io, print_str)
end

############################################################
@variables e f p q R C
r = DissipatorElement([e ~ R * f])
c = StaticStorageElement([e ~ C * q])

resistor = Component(r, "R1")
capacitor = Component(c, "C1")

j0 = Component(EqualEffort(), "j0")

b1 = Bond(resistor, j0)
b2 = Bond(j0, capacitor)

bg = BondGraph("NewBG", [resistor, capacitor, j0], [b1, b2])

# alternative construction
bg2 = BondGraph("NewBG", [b1, b2])

############################################################
# Graph representation
function graph(bg::BondGraph)
    bg_graph = MetaGraph(
        DiGraph();
        label_type=Symbol,
        vertex_data_type=Component,
        edge_data_type=Bond,
        graph_data=string(bg.name),
        # TODO add weight function and default weight
    )
    for comp in components(bg)
        bg_graph[comp.name] = comp
    end
    for bond in bg.bonds
        src, dst = vertices(bond)
        bg_graph[src.name, dst.name] = bond
    end
    bg_graph
end

g = graph(bg)
g[]

g[:C1]
g.vertex_properties
g.edge_data
g.vertex_labels

incidence_matrix(g)

############################################################################
using Plots, GraphRecipes
import GraphRecipes: graphplot

function graphplot(bg::BondGraph; kwargs...)
    g = graph(bg)
    graphplot(g;
        title = bg.name,
        names = g.vertex_labels,
        curves = false,
        nodeshape = :rect,
        kwargs...
    )
end

# graphplot(bg)

# TODO CONTINUE FROM HERE
