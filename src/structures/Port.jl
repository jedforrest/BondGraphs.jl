@connector function PortConnector(; name::Symbol)
    vars = @variables begin
        e(t), [connect = Effort]
        f(t), [connect = Flow]
    end
    return System(Equation[], t, vars, []; name)
end

# used for adjusting the sign of efforts/flows
# @enum PortWeight INBOUND=1 OUTBOUND=-1 DISCONNECTED=0

mutable struct Port
    name::Symbol
    parentname::Symbol
    model::System  # port connector
    connected::Bool
    # weight::PortWeight
    function Port(name, parentname)
        sys = PortConnector(; name)
        new(name, parentname, sys, false)
    end
end

name(p::Port) = p.name
parentname(p::Port) = p.parentname
function system(p::Port; namespaced = true)
    namespaced ? ModelingToolkitBase.renamespace(p.parentname, p.model) : p.model
end

is_connected(p::Port) = p.connected

is_source(p::Port) = !isnothing(p.inbound) && !p.inbound

portweight(p::Port) = isnothing(p.inbound) ? 0 : (p.inbound ? 1 : -1)

function Base.show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.parentname).$(port.name) $connection_state")
end
