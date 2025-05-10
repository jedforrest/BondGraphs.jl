import Base: show
using Symbolics

include("energypair.jl")

struct Port
    name::AbstractString
    variables::EnergyPair
    connected::Bool
end
function Port(name, e::Symbol=:e, f::Symbol=:f)
    Port(name, EnergyPair(e, f), false)
end
function show(io::IO, port::Port)
    is_connected = port.connected ? "⚫" : "⚪"
    print(io, "$(port.name) $is_connected")
end
variables(p::Port) = p.variables
is_connected(p::Port) = p.connected

port = Port("p1")
variables(port)

abstract type AbstractNode end

struct Component{C, N} <: AbstractNode where {N<:Integer}
    name::AbstractString
    ports::NTuple{N,Port}
    equations::Vector{Equation}
    function Component(class::Symbol, nports::Integer=1; name::AbstractString="")
        ports = ntuple(i -> Port("p$i"), nports)
        new{class,nports}(name, ports, Equation[])
    end
end
class(node::AbstractNode) = typeof(node).parameters[1]
name(node::AbstractNode) = node.name
numports(node::AbstractNode) = length(node.ports)
variables(node::AbstractNode) = variables.(node.ports)

function show(io::IO, comp::Component)
    print(io, "$(class(comp)){$(numports(comp))}")
end


c = Component(:C)
r = Component(:R, 2)

variables(c)
r.ports
numports(r)
variables(r)

struct Junction{C} <: AbstractNode
    name::AbstractString
    ports::Vector{Port}
    function Junction(class::Symbol; name::AbstractString="")
        new{class}(name, Port[])
    end
end
show(io::IO, junc::Junction) = print(io, "$(class(junc))")
Junction(:J0)


module Standard
# library of standard bond graph component definitions
# mostly one-port linear components

# Unlike the current approach, these should be functions
# that create instance of the structs themselves,
# rather than a dict or other description of a component

end
