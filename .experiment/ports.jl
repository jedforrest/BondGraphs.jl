import Base: show
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# mutable so that connection bool can be changed
mutable struct Port
    name::Symbol  # TODO index instead of name?
    # parent::Symbol
    effortflow::EffortFlowPair
    connected::Bool
    function Port(name, effortflow)
        new(name, effortflow, false)
    end
end
Port(e::Num, f::Num; name=:port) = Port(Symbol(name), EffortFlowPair(e, f))
Port(; name=:port) = Port(Symbol(name), EffortFlowPair())

effort(p::Port) = effort(p.effortflow)
flow(p::Port) = flow(p.effortflow)
power(p::Port) = effort(p) * flow(p)

name(p::Port) = p.name
parent(p::Port) = p.parent
is_connected(p::Port) = p.connected

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⚫" : "⚪"
    print(io, "$(name(port)) $connection_state")
end
