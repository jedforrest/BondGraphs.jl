"""
    Component{T<:AbstractElement}
    Component(type, name=type)
    Component(type, name=type; library=BondGraphs.DEFAULT_LIBRARY, <keyword arguments>)

Construct a Component of a defined (bondgraph) type ∈ {R, C, I, Se, Sf, TF, Ce, Re, SCe}.

Components have a `N` fixed ports when generated. This is usually determined by the bond
graph type. Other properties and equations of available components are defined in
`BondGraphs.DEFAULT_LIBRARY` (see  [`description`](@ref)).
"""
struct Component{T<:AbstractElement}
    element::T
    name::Symbol
    model::System
    efforts::Vector
    flows::Vector
    states::Vector
    ports::Vector{Port}
end

"""Generic Component constructor from a vector of equations"""
function Component(element::AbstractElement, eqs::Vector{Equation}, args...; name, kwargs...)
    sys = System(eqs, t; name=Symbol("$name.model"))
    Component(element, sys, args...; name, kwargs...)
end

"""Bond Element Component constructor from a MTK System"""
function Component(element::AbstractElement, sys::System, efforts, flows, states=Num[]; name=nameof(sys), kwargs...)
    check_num_ports(element, efforts, flows, states)

    # make new port subsystems that expose the effort and flow variables
    ports = Port[]
    portsystems = System[]
    port_connection_eqs = Equation[]
    for i in 1:numports(element)
        # create N port subsystems and extend the user-given MTK System
        port = Port(Symbol("port_$i"), name)
        portsys = system(port, namespaced=false)

        # add effort/flow connections to newly added port variables
        # (assuming efforts and flows are in the correct order)
        port_conn_eq = [
            LocalScope(efforts[i]) ~ ParentScope(portsys.e),
            LocalScope(flows[i]) ~ ParentScope(portsys.f)
        ]
        # append!(port_connection_eqs, port_conn_eq)

        # extend the port subsystem instead
        # portsys = extend(portsys, System(port_conn_eq, t; name))

        sys = compose(sys, portsys)
        sys = extend(sys, System(port_conn_eq, t; name))

        push!(ports, port)
        push!(portsystems, portsys)
    end
    # sys = extend(sys, System(port_connection_eqs, t; name))
    # sys = compose(sys, portsystems)

    # kwargs are used to set default parameter values
    for (key, val) in kwargs
        if hasproperty(sys, key)
            ModelingToolkitBase.get_initial_conditions(sys)[getproperty(sys, key)] = val
        end
    end

    Component(element, name, sys, efforts, flows, states, ports)
end

"""Non parameteric junction constructor"""
function Component(element::NonParametricJunction; name::Symbol)
    sys = System(Equation[], t; name)
    Component(element, name, sys, Num[], Num[], Num[], Port[])
end


############################################################
function check_num_ports(element::BondElement, efforts, flows, states)
    nports = numports(element)
    nports >= 1 || error("Must have at least 1 port ($nports)")
    if element isa StorageElement
        (nports == length(efforts) == length(flows) == length(states)) ||
            error("Number of efforts, flows, and states don't match ($nports)")
    else
        (nports == length(efforts) == length(flows)) ||
            error("Number of efforts and flows don't match ($nports)")
    end
    return nothing
end
function check_num_ports(element::SourceElement, efforts, flows, states)
    nports = numports(element)
    (nports == length(efforts) == length(flows) == 1) ||
        error("Must have exactly 1 port ($nports)")
    return nothing
end
function check_num_ports(element::ParametricJunction, efforts, flows, states)
    nports = numports(element)
    nports >= 2 || error("Must have at least 2 ports ($nports)")
    (nports == length(efforts) == length(flows)) ||
        error("Number of efforts and flows don't match ($nports)")
    return nothing
end
function check_num_ports(::NonParametricJunction, efforts, flows, states)
    # nothing to check
    return nothing
end

############################################################
# PROPERTIES
elementtype(::Component{T}) where {T} = T
element(comp::Component) = comp.element

name(comp::Component) = comp.name

symbol(comp::Component) = symbol(comp.element)

system(comp::Component) = comp.model
# System construction for junction components
# Since the number of ports can change, these systems are composed each time they are called
function system(junc0::Component{ZeroJunction})
    isempty(junc0.ports) && return junc0.model

    _ports = ports(junc0)
    es = effort.(_ports)
    fs = flow.(_ports)

    # e_in = -e_out
    # Σ f_in = Σ f_out
    eqs = [
        [es[1] ~ e for e in es[2:end]];
        sum(fs; init=0) ~ 0
    ]
    System(eqs, t; systems=system.(_ports), name=name(junc0))
end
function system(junc1::Component{OneJunction})
    isempty(junc1.ports) && return junc1.model

    _ports = ports(junc1)
    W = portweight.(_ports)
    es = W .* effort.(_ports)
    fs = W .* flow.(_ports)

    # Σ e_in = Σ e_out
    # f_in = -f_out
    eqs = [
        sum(es; init=0) ~ 0;
        [fs[1] ~ f for f in fs[2:end]]
    ]
    System(eqs, t; systems=system.(_ports), name=name(junc1))
end

variables(comp::Component) = ModelingToolkitBase.unknowns_toplevel(system(comp))
parameters(comp::Component) = ModelingToolkitBase.parameters(system(comp))

efforts(comp::Component) = comp.efforts
flows(comp::Component) = comp.flows
states(comp::Component) = comp.states

function constitutive_relations(comp::Component)
    sys = system(comp)
    top_level = [
        ModelingToolkitBase.unknowns_toplevel(sys);
        ModelingToolkitBase.parameters_toplevel(sys)
    ]
    eqs = ModelingToolkitBase.equations_toplevel(sys)
    return filter(eqs) do eq
        # check non-derivative variables only
        eq_vars = filter(!Symbolics.is_derivative, ModelingToolkitBase.get_variables(eq))
        all(any(isequal(v), top_level) for v in eq_vars)
    end
end

ports(comp::Component) = comp.ports

##############################
# BASE FUNCTIONS

function Base.show(io::IO, comp::Component)
    print(io, "$(numports(comp))-port $(symbol(comp))-Component $(name(comp))")
end

function Base.show(io::IO, comp::Component{<:NonParametricJunction})
    print(io, "$(symbol(comp))-Junction $(name(comp))")
end

# Easier referencing systems using a.b notation
function Base.getproperty(comp::Component, name::Symbol)
    if isdefined(comp, name)
        return getfield(comp, name)
    else
        # get default value for variable/parameter if it exists
        sym = getproperty(comp.model, name; namespace=false)
        return get(defaults(comp.model), sym, nothing)
    end
end

function Base.setproperty!(comp::Component{<:BondElement}, name::Symbol, val)
    setproperty!(system(comp), name, val)
end

# index referencing for ports
Base.getindex(comp::Component, index::Int) = comp.ports[index]
function Base.getindex(comp::Component, key::Symbol)
    ports = comp.ports
    port_index = findfirst(x -> name(x) == key, ports)
    !isnothing(port_index) ? ports[port_index] : error("No such port: $key")
end

############################################################
numports(comp::Component) = length(ports(comp))

hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

function nextfreeport(comp::Component)
    freeports = filter(!is_connected, comp.ports)
    isempty(freeports) ? nothing : first(freeports)
end
function nextfreeport(junc::Component{<:NonParametricJunction})
    # 0- and 1- junctions have unlimited ports
    # so create a new port if none are free
    freeportindex = findfirst(!is_connected, junc.ports)
    if isnothing(freeportindex)
        index = length(junc.ports) + 1
        port = Port(Symbol("port_$index"), name(junc))
        push!(junc.ports, port) # add new port to junction
        return port
    else
        return junc.ports[freeportindex]
    end
end
