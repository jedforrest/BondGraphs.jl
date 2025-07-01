using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using StaticArrays

import Base: show

include("ontology.jl")
include("standardlibrary.jl")
using .Library
using .Library: PowerPort


@named rec = Library.Reaction()

function is_powerport(model::ModelingToolkit.AbstractSystem)
    ModelingToolkit.isconnector(model) &&
        Set(getconnect.(unknowns(model))) == Set([Effort, Flow])
end

get_powerports(model::ModelingToolkit.AbstractSystem) = filter(is_powerport, ModelingToolkit.get_systems(model))

get_powerports(rec)

port = rec.systems[1]
is_powerport(port)

############################################################################

struct Port
    name::Symbol  # index instead of name?
    model::ModelingToolkit.AbstractSystem
    is_connected::Base.RefValue{Bool}
    function Port(model::ModelingToolkit.AbstractSystem)
        @assert is_powerport(model)
        new(model.name, model, Ref(false))
    end
end

# assuming vars are labelled :e and :f
effort(p::Port) = p.model.e[1]
flow(p::Port) = p.model.f[1]
power(p::Port) = effort(p) * flow(p)

name(p::Port) = p.name
model(p::Port) = p.model
is_connected(p::Port) = p.is_connected[]

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(name(port)) $connection_state")
end

p = Port(port)

############################################################################
# TODO CONTINUE FROM HERE
# - ports and MTK models have been created, now create Elements for the bond graph
# - will need to somehow add models in a way consistent with the BondGraphElement Hierarchy
#   - possibly using dispatched constructors

struct Element{T<:BondGraphElement}
    name::Symbol
    ports::Vector{Port}
    model::Union{ModelingToolkit.Model,Nothing}  # Optional at init
end

function Element{T}(model=nothing; name) where {T<:BondGraphElement}
    ports = isnothing(model) ? Port[] : Port.(get_powerports(model))
    Element{T}(Symbol(name), ports, model)
end

name(node::Element) = node.name
numports(node::Element) = length(node.ports)
model(node::Element) = node.model
class(node::Element) = typeof(node).parameters[1]

hasmodel(node::Element) = !isnothing(model(node))


# Used when displaying in a graph.
glyph(::Element) = :X
glyph(::Element{BondGraphElement}) = :E

glyph(::Element{StaticStorageElement}) = :C
glyph(::Element{DynamicStorageElement}) = :I
glyph(::Element{DissipatorElement}) = :R
# glyph(::Component{EffortSource}) = :Se
# glyph(::Component{FlowSource}) = :Sf
# glyph(::Component{Transformer}) = :TF
# glyph(::Component{Gyrator}) = :GY

# glyph(::Component{JunctionStructure}) = :J
# glyph(::Component{EqualEffort}) = Symbol(0)
# glyph(::Component{EqualFlow}) = Symbol(1)


show(io::IO, comp::Element) = print(io, "$(glyph(comp))::$(name(comp))")
show(io::IO, comp::Element{<:JunctionStructure}) = print(io, "$(glyph(comp))")


el = Element{StaticStorageElement}(name="foo")
class(el)
numports(el)
hasmodel(el)
glyph(el)

############################################################################

# hasfreeport(comp::Component) = any(!is_connected, ports(comp))
# hasfreeport(::Component{<:NonParametricJunction}) = true

# nextfreeport(comp::Component) = first(filter(!is_connected, ports(comp)))
# function nextfreeport(junc::Component{<:NonParametricJunction})
#     # FIXME should only create ports if none are free
#     # 0- and 1- junctions have unlimited ports
#     # so create a new port if trying to connect
#     index = numports(junc) + 1
#     port = Port(name="p$index")
#     push!(junc.ports, port)
#     port
# end


############################

# TODO determine whether a defined MTK.AbstractSystem fits the component definition
# e.g. a storage or dissipative system is correctly described as such
