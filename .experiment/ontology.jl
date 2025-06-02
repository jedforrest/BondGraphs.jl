# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

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

# Julia Type trees
# using GraphRecipes, Plots
# default(size=(1000, 1000))
# plot(BondGraphVertexClass, method=:tree, fontsize=10, nodeshape=:ellipse)

############################################################################

# abstract type AbstractNode end  # with the new ontology, this may not be needed?

# FIXME CONTINUE FROM HERE
# Challenge: how to test that MTK models are chosen correctly
# i.e. that a storage component is really a storage component

struct Component{T<:BondGraphVertexClass}
    class::T
    name::Symbol
    model::ModelingToolkit.AbstractSystem
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

############################################################################

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


############################
# TODO determine whether a defined MTK.AbstractSystem fits the component definition
# e.g. a storage or dissipative system is correctly described as such
using ModelingToolkitStandardLibrary

using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant

typeof(__Resistor__)
ModelingToolkit.Model

@named c = Capacitor()
hierarchy(c)
typeof(Capacitor)
typeof(c)

Capacitor.structure
c
equations(expand_connections(c))
equations(c)
