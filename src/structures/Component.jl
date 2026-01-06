struct Port
    name::Symbol
    sys::System  # port subsys
    parentname::Symbol
    connected::Ref{Bool}
    function Port(name, sys, parentname)
        new(name, sys, parentname, Ref(false))
    end
end
Port(sys::System, parentname::Symbol) = Port(nameof(sys), sys, parentname)

name(p::Port) = p.name
function system(p::Port; namespaced = true)
    namespaced ? ModelingToolkit.renamespace(p.parentname, p.sys) : p.sys
end
parent(p::Port) = p.parentname

is_connected(p::Port) = p.connected[]
connect!(p::Port) = p.connected[] = true

function Base.show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.parentname).$(port.name) $connection_state")
end

############################################################
abstract type AbstractNode end

"""
    Component{N} <: AbstractNode
    Component(type, name=type)
    Component(type, name=type; library=BondGraphs.DEFAULT_LIBRARY, <keyword arguments>)

Construct a Component of a defined (bondgraph) type ∈ {R, C, I, Se, Sf, TF, Ce, Re, SCe}.

Components have a `N` fixed ports when generated. This is usually determined by the bond
graph type. Other properties and equations of available components are defined in
`BondGraphs.DEFAULT_LIBRARY` (see  [`description`](@ref)).
"""
struct Component{T<:BondElement} <: AbstractNode
    element::T
    name::Symbol
    sys::System
    ports::Vector{Port}
end

function Component(element::BondGraphVertex; name::Symbol)
    # create MTK System from element definition
    sys = system(element; name)
    powerports = filter(ModelingToolkit.isconnector, ModelingToolkit.get_systems(sys))
    bg_ports = isempty(powerports) ? Port[] : Port.(powerports, name)
    Component(element, name, sys, bg_ports)
end

############################################################

struct Junction{T<:JunctionStructure} <: AbstractNode
    element::T
    name::Symbol
    sys::System
    ports::Vector{Port}
end

"""
    EqualEffort <: Junction

Efforts are all equal, flows sum to zero (0-junction). Has an unlimited number of ports.
"""

"""
    EqualFlow <: Junction

Flows are all equal, efforts sum to zero (1-junction). Has an unlimited number of ports.
"""

############################################################
# PROPERTIES
elementtype(::Component{V}) where {V} = V

variables(comp::Component) = ModelingToolkit.get_unknowns(system(comp))  # FIXME should be toplevel only
parameters(comp::Component) = ModelingToolkit.parameters(system(comp))

efforts(comp::Component) = getefforts(system(comp))
flows(comp::Component) = getflows(system(comp))

constitutive_relations(comp::Component) = equations(system(comp))  # FIXME should be toplevel only

############################################################
hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

function nextfreeport(comp::Component)
    freeports = filter(!is_connected, comp.ports)
    isempty(freeports) ? nothing : first(freeports)
end
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    portsys = PowerPort(; name = Symbol("port_$index"))
    port = Port(portsys, junc.name)
    push!(junc.ports, port) # add new port to junction
    port
end

##############################
system(comp::Component) = comp.sys

# System construction for junction components
# Since the number of ports can change, these systems are composed each time they are called
function system(junc::Component{<:NonParametricJunction})
    # create base system + port subsystems
    basesys = junc.sys
    isempty(junc.ports) && return basesys

    # base system without internal connection equations
    portsys = system.(junc.ports, namespaced = false)
    sys = compose(basesys, portsys)

    # create effort and flow conservation laws
    effort_eqs = junc.element.effort_fn(getefforts(sys))
    flow_eqs = junc.element.flow_fn(getflows(sys))
    extend(System([effort_eqs; flow_eqs], t; name = nameof(sys)), sys)
end


# Type
# type(n::AbstractNode) = n.type
# type(j::Junction) = typeof(j)
# type(::SourceSensor) = "SS"

# # Name
# name(n::AbstractNode) = n.name
# name(n::Junction) = vertex(n) == 0 ? n.name : "$(n.name)_$(vertex(n))"

# # Ports
# ports(n::AbstractNode) = n.ports
# numports(n::AbstractNode) = length(ports(n))
# portlabels(n::AbstractNode) = collect(keys(ports(n)))

# isconnected(n::AbstractNode, label) = ports(n)[label] != 0 # '!=0' needed for junctions
# updateport!(n::AbstractNode, label) = ports(n)[label] = !ports(n)[label]
# updateport!(::Junction, ::Int) = nothing # override

# port_info(n::AbstractNode) = (n, nextfreeport(n))
# port_info(t::Tuple{AbstractNode, Any}) = t

# # ports renamed as ports to make purpose clearer
# @deprecate freeports(n::AbstractNode) ports(n::AbstractNode)

# # Weights
# weights(j::Junction) = j.ports # deprecated
# set_weight!(j::Junction, idx::Int, w::Int) = ports(j)[idx] = w

# nextfreeport(n::AbstractNode) = findfirst(!, ports(n)) # first 'not' connected port
# function nextfreeport(j::Junction)
#     index = findfirst(==(0), ports(j))
#     if isnothing(index)
#         push!(j.ports, 0) # add new empty port
#         return numports(j)
#     else
#         return index
#     end
# end

# # Vertex
# vertex(n::AbstractNode) = n.vertex[]
# set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# # Parameters
# parameters(::AbstractNode) = Dict()
# parameters(n::Component) = n.variables[:parameters]

# # Globals
# globals(::AbstractNode) = Dict()
# globals(n::Component) = n.variables[:globals]

# # State variables
# states(::AbstractNode) = Dict()
# states(n::Component) = n.variables[:states]

# # Control variables
# controls(::AbstractNode) = Dict()
# controls(n::Component) = n.variables[:controls]

# # Equations
# equations(::AbstractNode) = Equation[]
# equations(n::Component) = n.equations

# # Variables
# all_variables(::AbstractNode) = ()
# function all_variables(n::Component)
#     # use getfield here so that this can be used by getproperty
#     merge(values(getfield(n, :variables))...)
# end

##############################
# BASE FUNCTIONS
# This definition will need to expand when equations etc. are added
Base.(==)(n1::AbstractNode, n2::AbstractNode) = type(n1) == type(n2) && name(n1) == name(n2)

# Base.show(io::IO, node::AbstractNode) = print(io, "$(type(node)):$(name(node))")
# Base.show(io::IO, node::Junction) = print(io, name(node))

function Base.show(io::IO, comp::Component{<:BondElement})
    print(io, "$(glyph(comp.element))::$(comp.name)")
end
Base.show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.element))")



# Easier referencing systems using a.b notation
# TODO: rearrange so that getfield() is checked first with isdefined()
# function getproperty(n::Component, sym::Symbol)
#     p, = @parameters $sym
#     _, x = @variables t, $sym(t)
#     all_vars = all_variables(n)

#     if p in keys(all_vars)
#         return all_vars[p]
#     elseif x in keys(all_vars)
#         return all_vars[x]
#     else
#         getfield(n, sym)
#     end
# end

# # TODO: this can overrite global variables unintentionally by creating a new local copy
# function setproperty!(n::Component, sym::Symbol, val)
#     p, = @parameters $sym
#     _, x = @variables t, $sym(t)

#     for (_, vars) in getfield(n, :variables)
#         if p in keys(vars)
#             return vars[p] = val
#         elseif x in keys(vars)
#             return vars[x] = val
#         end
#     end
#     setfield!(n, sym, val)
# end
