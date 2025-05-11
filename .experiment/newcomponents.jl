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
effort(p::Port) = effort(p.variables)
flow(p::Port) = flow(p.variables)
is_connected(p::Port) = p.connected

port = Port("p1")
variables(port)
effort(port)
flow(port)

abstract type AbstractNode end

struct Component{C, N} <: AbstractNode where {N<:Integer}
    name::AbstractString
    ports::NTuple{N,Port}
    equations::Vector{Equation}
    function Component(class::Symbol, nports::Integer=1; name="", equations=Equation[])
        ports = ntuple(i -> Port("p$i"), nports)
        new{class,nports}(name, ports, equations)
    end
end
class(node::AbstractNode) = typeof(node).parameters[1]
name(node::AbstractNode) = node.name
numports(node::AbstractNode) = length(node.ports)
variables(node::AbstractNode) = variables.(node.ports)
# for now, equation and port variables are different
# it will be easier to connect variables using MTK.jl connectors (TODO)
equations(node::AbstractNode) = node.equations

function show(io::IO, comp::Component)
    print_str = "$(class(comp)){$(numports(comp))}"
    print_str *= isempty(name(comp)) ? "" : " [$(name(comp))]"
    print(io, print_str)
end

struct Junction{C} <: AbstractNode
    name::AbstractString
    ports::Vector{Port}
    function Junction(class::Union{Symbol,Integer}; name::AbstractString="")
        new{class}(name, Port[])
    end
end
show(io::IO, junc::Junction) = print(io, "$(class(junc)){$(numports(junc))}")

function add_port!(junc::Junction)
    index = numports(junc) + 1
    newport = Port("p$index", Symbol("e_$index"), Symbol("f_$index"))
    push!(junc.ports, newport)
end

function equations(junc::Junction{0})
    efforts, flows = effort.(junc.ports), flow.(junc.ports)
    effort_eqns = efforts[1] .~ efforts[2:end]
    flow_eqns = sum(flows) ~ 0
    [effort_eqns; flow_eqns]
end
function equations(junc::Junction{1})
    efforts, flows = effort.(junc.ports), flow.(junc.ports)
    effort_eqns = sum(efforts) ~ 0
    flow_eqns = flows[1] .~ flows[2:end]
    [effort_eqns; flow_eqns]
end

j0 = Junction(0)
j1 = Junction(1)

for i in 1:5
    add_port!(j0)
    add_port!(j1)
end
equations(j0)
equations(j1)


# module Standard
# library of standard bond graph component definitions
# mostly one-port linear components

# Unlike the current approach, these should be functions
# that create instance of the structs themselves,
# rather than a dict or other description of a component

# TODO construct using MTK tools

using Symbolics

@variables t e(t) f(t) q(t) p(t)
@variables C R

D = Differential(t)

function capacitor(name="cap")
    equations = [q ~ C * e, D(q) ~ f]
    return Component(:C, 1; name, equations)
end

function resistor(name="res")
    equations = [e ~ R * f]
    return Component(:R, 1; name, equations)
end

# export capacitor, resistor

# end

# using .Standard

# c = Standard.capacitor()
c = capacitor()
r = resistor("r_comp")

variables(c)
r.ports
numports(r)
variables(r)
