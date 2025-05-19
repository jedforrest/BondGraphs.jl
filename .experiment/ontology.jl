# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit

include("energypair.jl")
include("ports.jl")

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

abstract type BondGraphVertexClass end

abstract type BondGraphElement <: BondGraphVertexClass end
abstract type JunctionStructure <: BondGraphVertexClass end

abstract type StorageElement <: BondGraphElement end
abstract type DissipatorElement <: BondGraphElement end
abstract type SourceElement <: BondGraphElement end

abstract type StaticStorageElement <: StorageElement end
abstract type DynamicStorageElement <: StorageElement end

abstract type EffortSource <: SourceElement end
abstract type FlowSource <: SourceElement end

abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

abstract type Transformer <: ParametricJunction end
abstract type Gyrator <: ParametricJunction end

abstract type EqualEffort <: NonParametricJunction end
abstract type EqualFlow <: NonParametricJunction end

############################################################################

abstract type AbstractNode end  # with the new ontology, this may not be needed?

struct Component{T<:BondGraphVertexClass} <: AbstractNode
    name::Symbol
    variables::Vector{Num}
    parameters::Vector{Num}
    equations::Vector{Equation}
    ports::Vector{Port}
end
function Component{T}(; name, variables=Num[], parameters=Num[], equations=Equation[], ports=Port[]) where {T<:BondGraphVertexClass}
    # empty 'generic' component
    Component{T}(name, variables, parameters, equations, ports)
end
name(node::Component) = node.name
variables(node::Component) = node.variables
parameters(node::Component) = node.parameters
equations(node::Component) = node.equations
ports(node::Component) = node.ports

class(node::Component) = typeof(node).parameters[1]
numports(node::Component) = length(ports(node))
efforts(node::Component) = effort.(ports(node))
flows(node::Component) = flow.(ports(node))

c = Component{BondGraphVertexClass}(name=:test, ports=[Port()])
class(c)
numports(c)
efforts(c)
flows(c)

# or icon; can't think of a better word. Used when displaying in a graph.
glyph(::Component) = :X
glyph(::Component{BondGraphElement}) = :E
glyph(::Component{JunctionStructure}) = :J

glyph(::Component{StaticStorageElement}) = :C
glyph(::Component{DynamicStorageElement}) = :I
glyph(::Component{DissipatorElement}) = :R
glyph(::Component{EffortSource}) = :Se
glyph(::Component{FlowSource}) = :Sf

glyph(::Component{Transformer}) = :TF
glyph(::Component{Gyrator}) = :GY
glyph(::Component{EqualEffort}) = :𝟘
glyph(::Component{EqualFlow}) = :𝟙

glyph(c)

show(io::IO, comp::Component) = print(io, "$(glyph(comp))::$(name(comp))")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp))")

############################################################################

hasfreeport(comp::Component) = any(!is_connected, ports(comp))
hasfreeport(::Component{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, ports(comp)))
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = numports(junc) + 1
    port = Port(name="p$index")
    push!(junc.ports, port)
    port
end

struct Bond
    src::Port
    dst::Port
    function Bond(src::Port, dst::Port)
        src.connected = true
        dst.connected = true
        new(src, dst)
    end
end
function Bond(src_comp::Component, dst_comp::Component)
    hasfreeport(src_comp) || error("$src_comp has no free ports")
    hasfreeport(dst_comp) || error("$dst_comp has no free ports")
    Bond(nextfreeport(src_comp), nextfreeport(dst_comp))
end

# show(io::IO, b::Bond) = print(io, "Bond()")  # FIXME

# is_connected.(capacitor.ports)

# hasfreeport(capacitor)
# hasfreeport(kvl)

# nextfreeport(capacitor)
# nextfreeport(kvl)

############################################################################
import Base: size

# New bond graph structure
mutable struct NewBondGraph
    name::Symbol
    graph::MetaGraph
    function NewBondGraph(graph::AbstractGraph; name::Symbol)
        # creating MetaGraph in a constructor keeps it type stable
        metagraph = MetaGraph(
            graph;  # underlying graph structure
            label_type=Symbol,  # node name
            vertex_data_type=AbstractNode,  # node type
            edge_data_type=Bond,  # bond
            graph_data=name,  # tag for the whole graph
        )
        return new(name, metagraph)
    end
end
NewBondGraph(; name=:NewBG) = NewBondGraph(DiGraph(); name=Symbol(name))

name(bg::NewBondGraph) = bg.name
size(bg::NewBondGraph) = (nv(bg.graph), ne(bg.graph))

show(io::IO, bg::NewBondGraph) = print(io, "$(name(bg)) BondGraph$(size(bg))")

bg = NewBondGraph()

############################################################################

function add_node!(bg::NewBondGraph, comp::Component)
    bg.graph[comp.name] = comp
end

function connect!(bg::NewBondGraph, src_comp::Component, dst_comp::Component)
    # TODO assuming single ports, change to allow selecting a specific port
    bg.graph[src_comp.name, dst_comp.name] = Bond(src_comp, dst_comp)
end

############################################################################
############################################################################

@variables t
D = Differential(t)

function Capacitor(; name=:capacitor)
    variables = @variables e(t), f(t), q(t)
    parameters = @parameters C
    equations = [q ~ C * e, D(q) ~ f]
    ports = [Port(e, f; name="pC")]
    return Component{StaticStorageElement}(; name, variables, parameters, equations, ports)
end

function Inductor(; name=:inductor)
    variables = @variables e(t), f(t), p(t)
    parameters = @parameters L
    equations = [p ~ L * f, D(p) ~ e]
    ports = [Port(e, f; name="pI")]
    return Component{DynamicStorageElement}(; name, variables, parameters, equations, ports)
end

function Resistor(; name=:resistor)
    variables = @variables e(t), f(t)
    parameters = @parameters R
    equations = [e ~ R * f]
    ports = [Port(e, f; name="pR")]
    return Component{DissipatorElement}(; name, variables, parameters, equations, ports)
end

function ZeroJunction(; name=:zero)
    # variables = @variables e(t), f(t)
    # ports = [Port(e, f)]
    return Component{EqualEffort}(; name)
end

############################################################################

@named rc_model = NewBondGraph()

capacitor = Capacitor()
resistor = Resistor()
kvl = ZeroJunction()

add_node!(rc_model, capacitor)
add_node!(rc_model, resistor)
add_node!(rc_model, kvl)

connect!(rc_model, capacitor, kvl)
connect!(rc_model, kvl, resistor)

rc_model

incidence_matrix(rc_model.graph)

# TODO CONTINUE
using Plots
plot(rc_model.graph)
