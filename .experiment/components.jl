using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using StaticArrays

import Base: show

# include("ontology.jl")
include("standardlibrary.jl")
using .Library
# using .Library: PowerPort


# @named rec = Library.Reaction()

function is_powerport(model::ModelingToolkit.AbstractSystem)
    ModelingToolkit.isconnector(model) &&
        Set(getconnect.(unknowns(model))) == Set([Effort, Flow])
end

get_powerports(model::ModelingToolkit.AbstractSystem) = filter(is_powerport, ModelingToolkit.get_systems(model))


############################################################################

const Effort = ModelingToolkit.Equality

@connector PowerPort begin
    e(t) = 0., [connect = Effort]
    f(t) = 0., [connect = Flow]
end

struct Port
    name::Symbol
    is_connected::Base.RefValue{Bool}
    Port(name) = new(Symbol(name), Ref(false))
end
Port(i::Integer) = Port("p$i")

# assuming vars are labelled :e and :f
# effort(p::Port) = p.sys.e
# flow(p::Port) = p.sys.f
# power(p::Port) = effort(p) * flow(p)

name(p::Port) = p.name
# model(p::Port) = p.sys
is_connected(p::Port) = p.is_connected[]

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(name(port)) $connection_state")
end

p = Port(1)

############################################################################
# TODO CONTINUE FROM HERE
# - ports and MTK models have been created, now create Elements for the bond graph
# - will need to somehow add models in a way consistent with the BondGraphElement Hierarchy
#   - possibly using dispatched constructors

# @component function ConstitutiveRelations(eqs, es, fs)
#     length(es) == length(fs) || throw(error("Number of efforts and flows must match"))
#     ODESystem(eqs, t, vars)
# end

############################################################################

# TODO CONTINUE FROM HERE
# likely need a CR struct so that we know what the effort and flows are

# mutable struct ConstitutiveRelation
#     eqs::Vector{Equation}
#     es::Vector{Num}
#     fs::Vector{Num}
#     function ConstitutiveRelation(eqs, es, fs)
#         # eqn_vars = Set(get_variables(eqn))
#         # for var in [es; fs]
#         #     var in eqn_vars || throw(error("$var not found in equation"))
#         # end
#         es = setmetadata.(es, ModelingToolkit.VariableConnectType, Effort)
#         fs = setmetadata.(fs, ModelingToolkit.VariableConnectType, Flow)
#         new(eqs, es, fs)
#     end
# end

# Set(vcat(get_variables.([e~f, e~2*f]))...)

# show(io::IO, cr::ConstitutiveRelation) = print(io, cr.eqs)


############################################################################

struct Element{T<:BondElement}
    sys::ModelingToolkit.AbstractSystem
    ports::AbstractVector{Port}
end

function (BGE::Type{<:BondElement})(sys; nports=1)
    ports = MVector{nports}([Port(i) for i in 1:nports])
    Element{BGE}(sys, ports)
end

@named C = Capacitor()
StaticStorageElement(C)


function Element(cr::ConstitutiveRelation; name, nports=1)
    basesys = ODESystem
    Element{T}(sys, ports)
end

function Element{T}(model=nothing; name) where {T<:BondElement}
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
glyph(::Element{BondElement}) = :E

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


############################

# TODO determine whether a defined MTK.AbstractSystem fits the component definition
# e.g. a storage or dissipative system is correctly described as such
