# Elements define a particular type of bond graph component (e.g. an electrical capacitor).
# An element will have constitutive relations (equations), efforts/flows, and any states.
# A bond graph can have multiple instances of an element.
# See Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

############################################################
abstract type AbstractElement end
abstract type BondElement <: AbstractElement end

"""`R` component"""
struct DissipatorElement <: BondElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
end

#########################
abstract type StorageElement <: BondElement end

"""`C` component"""
struct StaticStorageElement <: StorageElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    states::Vector
    ports::Vector{Port}
end

"""`I` component"""
struct DynamicStorageElement <: StorageElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    states::Vector
    ports::Vector{Port}
end

############################################################
abstract type SourceElement <: BondElement end

"""`Se` component"""
struct EffortSource <: SourceElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
end

"""`Sf` component"""
struct FlowSource <: SourceElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
end

"""`SS` component"""
struct SourceSensor <: SourceElement
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
end

############################################################

function (TBE::Type{<:BondElement})(sys::System, efforts, flows, states=[]; name::Symbol=nameof(sys), kwargs...)
    check_num_ports(efforts, flows)

    # make new port subsystems that expose the effort and flow variables

    # create ports and add port connections
    ports = Port[]
    portsystems = System[]
    port_connection_eqs = Equation[]
    for (i, (e, f)) in enumerate(zip(efforts, flows))
        # create N port subsystems and extend the user-given MTK System
        port = Port(Symbol("port_", i), name)
        portsys = system(port, namespaced=false)

        # add effort/flow connections to newly added port variables
        # (assuming efforts and flows are in the correct order)
        port_conn_eq = [
            ParentScope(e) ~ ParentScope(portsys.e),
            ParentScope(f) ~ ParentScope(portsys.f)
        ]

        push!(ports, port)
        push!(portsystems, portsys)
        append!(port_connection_eqs, port_conn_eq)
    end
    sys = compose(sys, portsystems)
    sys = compose(sys, System(port_connection_eqs, t; name))

    # kwargs are used to set default parameter values
    for (key, val) in kwargs
        setproperty!(sys, key, val)
    end

    if hasfield(TBE, :states)
        return TBE(name, sys, efforts, flows, states, ports)
    else
        return TBE(name, sys, efforts, flows, ports)
    end
end

""" Generic BondGraph Element constructor from a vector of equations"""
function (TE::Type{<:AbstractElement})(eqs::Vector{Equation}, args...; kwargs...)
    sys = System(eqs, t; name=Symbol(TE))
    TE(sys, args...; kwargs...)
end

############################################################
abstract type JunctionStructure <: AbstractElement end
abstract type ParametricJunction <: JunctionStructure end

# Junctions are a power-conserving transformation between the effort/flows of each port

"""`TF` component"""
struct Transformer <: ParametricJunction
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
    # function Transformer(sys::System, efforts, flows)
    #     numports = length(efforts)
    #     numports >= 2 ||
    #         error("Transformer must have at least 2 ports ($numports ports given)")
    #     new(sys, efforts, flows)
    # end
end
"""`GY` component"""
struct Gyrator <: ParametricJunction
    name::Symbol
    sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
    # function Gyrator(sys::System, efforts, flows)
    #     numports = length(efforts)
    #     numports >= 2 || error("Gyrator must have at least 2 ports ($numports ports given)")
    #     new(sys, efforts, flows)
    # end
end

############################################################
abstract type NonParametricJunction <: JunctionStructure end

# Since the number of ports is unknown and can change,
# we instead store functions that map efforts and flows to equations
"""`0`-junction"""
struct EqualEffort <: NonParametricJunction
    name::Symbol
    sys::System
    relation::Function
    ports::Vector{Port}
    function EqualEffort(; name=:𝟎)
        sys = System(Equation[], t; name)
        relation_fn(es, fs) = [
            [es[1] ~ e for e in es[2:end]];
            sum(fs) ~ 0
        ]
        new(name, sys, relation_fn, Port[])
    end
end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction
    name::Symbol
    sys::System
    relation::Function
    ports::Vector{Port}
    function EqualFlow(; name=:𝟏)
        sys = System(Equation[], t; name)
        relation_fn(es, fs) = [
            sum(es) ~ 0;
            [fs[1] ~ f for f in fs[2:end]]
        ]
        new(name, sys, relation_fn, Port[])
    end
end

# function (TNPJ::Type{<:NonParametricJunction})(; name=Symbol(TNPJ))
#     sys = System(Equation[], t; name)
#     TNPJ(name, sys)
# end

############################################################


function check_num_ports(efforts, flows)
    numports = length(efforts)
    (numports == length(flows)) ||
        error("Number of efforts and flows don't match ($numports)")
    numports >= 1 || error("Must have at least 1 port ($numports)")
    return nothing
end

############################################################

name(elem::AbstractElement) = elem.name

efforts(elem::BondElement) = elem.efforts
efforts(junc::JunctionStructure) = effort.(ports(junc))
flows(elem::BondElement) = elem.flows
flows(junc::JunctionStructure) = flow.(ports(junc))

states(elem::StorageElement) = elem.states

equations(elem::AbstractElement) = equations(elem.sys)

ports(elem::AbstractElement) = elem.ports

# excludes equation relating to port connections
function constitutive_relations(elem::BondElement)
    ModelingToolkit.equations_toplevel(elem.sys)
end
function constitutive_relations(junc::JunctionStructure)
    junc.relation(efforts(junc), flows(junc))
end

numports(elem::AbstractElement) = length(elem.ports)
numports(::NonParametricJunction) = Inf

# Used when displaying in a graph.
# TODO these can just be included in the struct definitions above (kwdef)
glyph(::DissipatorElement) = :R
glyph(::StaticStorageElement) = :C
glyph(::DynamicStorageElement) = :I
glyph(::EffortSource) = :Se
glyph(::FlowSource) = :Sf
glyph(::Transformer) = :TF
glyph(::Gyrator) = :GY
glyph(::JunctionStructure) = :J
glyph(::EqualEffort) = :𝟎
glyph(::EqualFlow) = :𝟏

############################################################
# Overloading Base

function Base.show(io::IO, elem::T) where {T <: AbstractElement}
    print_str = "$(glyph(elem))::$(name(elem))"
    print(io, print_str)
end
Base.show(io::IO, ::T) where {T <: NonParametricJunction} = print(io, T)

# Easier referencing systems using a.b notation
function Base.getproperty(elem::AbstractElement, name::Symbol)
    if isdefined(elem, name)
        return getfield(elem, name)
    else
        # get default value for variable/parameter if it exists
        sym = getproperty(elem.sys, name; namespace=false)
        return get(defaults(elem.sys), sym, nothing)
    end
end

############################

# index referencing for ports
Base.getindex(elem::AbstractElement, index::Int) = elem.ports[index]

function Base.getindex(elem::AbstractElement, key::Symbol)
    ports = elem.ports
    port_index = findfirst(x -> x.name == key, ports)
    !isnothing(port_index) ? ports[port_index] : error("No such port: $key")
end

############################################################
hasfreeport(elem::AbstractElement) = any(!is_connected, elem.ports)
hasfreeport(::NonParametricJunction) = true

function nextfreeport(elem::AbstractElement)
    freeports = filter(!is_connected, elem.ports)
    isempty(freeports) ? nothing : first(freeports)
end
function nextfreeport(junc::NonParametricJunction)
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    port = Port(Symbol("port_$index"), junc.name)
    push!(junc.ports, port) # add new port to junction
    port
end
