using ModelingToolkit: get_connection_type, renamespace

const Effort = ModelingToolkit.Equality

function power_variables(;
        e = :e,
        f = :f,
        p = :p,
        q = :q)
    @variables begin
        $(Symbol(e))(t), [connect = Effort]
        $(Symbol(f))(t), [connect = Flow]
        $(Symbol(p))(t)
        $(Symbol(q))(t)
    end
end

@connector function PortConnector(; name, e=:e, f=:f)
    vars = @variables begin
        $(Symbol(e))(t), [connect = Effort]
        $(Symbol(f))(t), [connect = Flow]
    end
    return System(Equation[], t, vars, []; name)
end

############################################################

mutable struct Port
    name::Symbol
    parentname::Symbol
    sys::System  # port connector
    connected::Bool
    is_source::Bool  # true: source, false: destination (or disconnected)
    function Port(name, parentname; e=:e, f=:f)
        sys = PortConnector(; name, e, f)
        new(name, parentname, sys, false, false)
    end
end

name(p::Port) = p.name
parentname(p::Port) = p.parentname
function system(p::Port; namespaced = true)
    namespaced ? renamespace(p.parentname, p.sys) : p.sys
end

is_connected(p::Port) = p.connected

function port_weight(p::Port)
    if !p.connected
        return 0
    elseif p.is_source
        return -1  # source port
    else
        return 1  # destination port
    end
end

function effort(p::Port)
    e = filter(x -> get_connection_type(x) == Effort, unknowns(p.sys))[]
    renamespace(p.name, e)
end
function flow(p::Port)
    f = filter(x -> get_connection_type(x) == Flow, unknowns(p.sys))[]
    renamespace(p.name, f)
end

function Base.show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.parentname).$(port.name) $connection_state")
end
