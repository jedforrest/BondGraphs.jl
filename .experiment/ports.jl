import Base: show
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

const Effort = ModelingToolkit.Equality

@connector function MTKPort(e::Symbol=:e, f::Symbol=:f; name)
    vars = @variables begin
        $e(t), [connect = Effort]
        $f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, Num[]; name)
end

struct Port
    name::Symbol  # TODO index instead of name?
    effort::Num
    flow::Num
    model::ModelingToolkit.AbstractSystem  # may have to change this to MTK.Model
    connected::Base.RefValue{Bool}
    function Port(e::Symbol=:e, f::Symbol=:f; name=:port)
        model = MTKPort(e, f; name)
        effort = getproperty(model, e)
        flow = getproperty(model, f)
        new(name, effort, flow, model, Ref(false))
    end
end

effort(p::Port) = p.effort
flow(p::Port) = p.flow
power(p::Port) = effort(p) * flow(p)

name(p::Port) = p.name
model(p::Port) = p.model
is_connected(p::Port) = p.connected[]

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⚫" : "⚪"
    print(io, "$connection_state [$(port.effort), $(port.flow)]")
end


port = Port()


@named begin
    port_a = Port(:u, :v)
    port_b = Port(:V, :I)
end
