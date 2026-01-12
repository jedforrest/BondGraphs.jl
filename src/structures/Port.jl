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

struct Port
    name::Symbol
    parentname::Symbol
    sys::System  # port connector
    connected::Ref{Bool}
    function Port(name, parentname; e=:e, f=:f)
        sys = PortConnector(; name, e, f)
        new(name, parentname, sys, Ref(false))
    end
end

name(p::Port) = p.name
parentname(p::Port) = p.parentname
function system(p::Port; namespaced = true)
    namespaced ? renamespace(p.parentname, p.sys) : p.sys
end

is_connected(p::Port) = p.connected[]
connect!(p::Port) = p.connected[] = true
disconnect!(p::Port) = p.connected[] = false

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
