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

"""`R` component"""
struct DissipatorElement <: BondElement
    cr::Vector{Equation}
end

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

numports(C::BondElement) = length(C.cr)

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
    name::Symbol
    parent::Any
    connected::Bool
    Port(name, parent) = new(Symbol(name), parent, false)
end
is_connected(p::Port) = p.connected

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.name) $connection_state")
end

############################################################
# TODO Component and Junction structure could be combined (into e.g. "Element")
# if they end up being similar enough

struct Component{E<:BondElement}
    type::E
    name::Symbol
    ports::Vector{Port}
end
function Component(type::BondElement, name)
    newcomp = Component(type, Symbol(name), Port[])
    for i in 1:numports(type)
        push!(newcomp.ports, Port("p$i", newcomp))
    end
    newcomp
end

struct Junction{J<:JunctionStructure}
    type::J
    name::Symbol
    ports::Vector{Port}
end
function Junction(type::JunctionStructure, name)
    Junction(type, Symbol(name), Port[])
end

const ComponentOrJunction = Union{Component,Junction}

show(io::IO, comp::Component) = print(io, "$(glyph(comp.type))::$(comp.name)")
show(io::IO, comp::Junction) = print(io, "$(glyph(comp.type))")

hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Junction{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, comp.ports))
function nextfreeport(junc::Junction{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    port = Port("p$index", junc)
    push!(junc.ports, port)
    port
end

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
function Bond(src::ComponentOrJunction, dst::ComponentOrJunction)
    hasfreeport(src) || error("$src has no free ports")
    hasfreeport(dst) || error("$dst has no free ports")
    srcport = nextfreeport(src)
    dstport = nextfreeport(dst)
    Bond(srcport, dstport)
end

srcvertex(b::Bond) = b.src.parent
dstvertex(b::Bond) = b.dst.parent

show(io::IO, b::Bond) = print(io, "$(srcvertex(b)) ⇀ $(dstvertex(b))")

############################################################

struct BondGraph
    name::Symbol
    components::Vector{Component}
    junctions::Vector{Junction}
    bonds::Vector{Bond}
    function BondGraph(name, components=[], junctions=[], bonds=[])
        new(Symbol(name), components, junctions, bonds)
    end
end
# 3 argument constructor
function BondGraph(name, components_junctions=[], bonds=[])
    components = filter(x -> x isa Component, components_junctions)
    junctions = filter(x -> x isa Junction, components_junctions)
    BondGraph(Symbol(name), components, junctions, bonds)
end
vertices(bg::BondGraph) = [bg.components; bg.junctions]

# maybe not as the default "show" option (summary print instead)
function show(io::IO, bg::BondGraph)
    print_str = """
    BondGraph $(bg.name)
    ----------------------
    $(join(bg.bonds,'\n'))
    """
    print(io, print_str)
end

############################################################
@variables e f p q R C
r = DissipatorElement([e ~ R * f])
c = StaticStorageElement([e ~ C * q])

resistor = Component(r, "resistor")
resistor.ports

capacitor = Component(c, "capacitor")

j0 = Junction(EqualEffort(), "j0")
# push!(j0.ports, Port("in"))
# push!(j0.ports, Port("out"))

b1 = Bond(resistor, j0)
b2 = Bond(j0, capacitor)

bg = BondGraph("NewBG", [resistor, capacitor, j0], [b1, b2])
vertices(bg)

c_j = [resistor, capacitor, j0]

############################################################
# Graph representation
function graph(bg::BondGraph)
    # creating MetaGraph in a function keeps it type stable
    bg_graph = MetaGraph(
        DiGraph();
        label_type=Symbol,
        vertex_data_type=ComponentOrJunction,
        edge_data_type=Bond,
        graph_data=string(bg.name),
        # TODO add weight function and default weight
    )
    for vertex in vertices(bg)
        bg_graph[vertex.name] = vertex
    end
    for bond in bg.bonds
        srcname = srcvertex(bond).name
        dstname = dstvertex(bond).name
        bg_graph[srcname, dstname] = bond
    end
    bg_graph
end

g = graph(bg)
g[]

g[:capacitor]
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

graphplot(bg)

# TODO CONTINUE FROM HERE
