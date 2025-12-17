# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
import Base: show, size
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations

############################################################
# abstract type AbstractPowerVariableType end

# struct EffortVar <: AbstractPowerVariableType end
# struct FlowVar <: AbstractPowerVariableType end
# struct MomentumVar <: AbstractPowerVariableType end
# struct DisplacementVar <: AbstractPowerVariableType end

############################################################
abstract type BondGraphVertex end

abstract type BondElement <: BondGraphVertex end
abstract type StorageElement <: BondElement end
abstract type SourceElement <: BondElement end

const Effort = ModelingToolkit.Equality

@connector PowerPort begin
    e(t), [connect = Effort]
    f(t), [connect = Flow]
end

############################################################
getefforts(sys::ODESystem) = [e for e in unknowns(sys) if getconnect(e) == Effort]
getflows(sys::ODESystem) = [f for f in unknowns(sys) if getconnect(f) == Flow]
getstates(sys::ODESystem) = [x for x in unknowns(sys) if !hasconnect(x)]

function check_num_ports(efforts, flows)
    numports = length(efforts)
    numports >= 1 || error("Must have at least 1 port ($numports)")
    (numports == length(flows)) || error("Number of efforts and flows don't match ($numports)")
    return nothing
end

############################################################
# TODO instead of storing a model or sys, this can store the symbolic equations?
# and label the efforts and flows (and state vars) explicitly
# NOTE: SDESystems do not yet work with @mtkmodel - therefore must use Non-DSL approach

"""`R` component"""
struct DissipatorElement <: BondElement
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    function DissipatorElement(eqs, efforts, flows)
        check_num_ports(efforts, flows)
        new(eqs, efforts, flows)
    end
end
function DissipatorElement(sys::ODESystem)
    DissipatorElement(equations(sys), getefforts(sys), getflows(sys))
end


"""`C` component"""
struct StaticStorageElement <: StorageElement
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    states::Vector
    function StaticStorageElement(eqs, efforts, flows, states)
        check_num_ports(efforts, flows)
        new(eqs, efforts, flows, states)
    end
end
function StaticStorageElement(sys::ODESystem)
    StaticStorageElement(equations(sys), getefforts(sys), getflows(sys), getstates(sys))
end


"""`I` component"""
struct DynamicStorageElement <: StorageElement
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    states::Vector
    function DynamicStorageElement(eqs, efforts, flows, states)
        check_num_ports(efforts, flows)
        new(eqs, efforts, flows, states)
    end
end
function DynamicStorageElement(sys::ODESystem)
    DynamicStorageElement(equations(sys), getefforts(sys), getflows(sys), getstates(sys))
end

############################################################
# TODO update structs with eqs, vars, etc.
"""`Se` component"""
struct EffortSource <: SourceElement
    eqs::Vector{Equation}
    efforts::Vector
    function EffortSource(eqs, efforts)
        numports = length(efforts)
        numports == 1 || error("Must have exactly 1 port ($numports)")
        new(eqs, efforts)
    end
end

"""`Sf` component"""
struct FlowSource <: SourceElement
    eqs::Vector{Equation}
    flows::Vector
    function FlowSource(eqs, flows)
        numports = length(flows)
        numports == 1 || error("Must have exactly 1 port ($numports)")
        new(eqs, flows)
    end
end

"""`SS` component"""
struct SourceSensor <: SourceElement
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    function SourceSensor(eqs, efforts, flows)
        numports = length(efforts)
        (numports == length(flows) == 1) || error("Must have exactly 1 port ($numports)")
        new(eqs, efforts, flows)
    end
end

############################################################
abstract type JunctionStructure <: BondGraphVertex end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

"""`TF` component"""
struct Transformer <: ParametricJunction
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    function Transformer(eqs, efforts, flows)
        numports = length(efforts)
        numports >= 2 || error("Transformer must have at least 2 ports ($numports ports given)")
        new(eqs, efforts, flows)
    end
end
"""`GY` component"""
struct Gyrator <: ParametricJunction
    eqs::Vector{Equation}
    efforts::Vector
    flows::Vector
    function Gyrator(eqs, efforts, flows)
        numports = length(efforts)
        numports >= 2 || error("Gyrator must have at least 2 ports ($numports ports given)")
        new(eqs, efforts, flows)
    end
end

############################################################

# Since the number of ports is unknown and can change,
# we instead store functions that map efforts/flows to equations
"""`0`-junction"""
struct EqualEffort <: NonParametricJunction
    effort_fn::Function
    flow_fn::Function
    function EqualEffort()
        effort_fn(es) = [es[1] ~ e for e in es[2:end]]
        flow_fn(fs) = sum(fs) ~ 0
        new(effort_fn, flow_fn)
    end
end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction
    effort_fn::Function
    flow_fn::Function
    function EqualFlow()
        effort_fn(es) = sum(es) ~ 0
        flow_fn(fs) = [fs[1] ~ f for f in fs[2:end]]
        new(effort_fn, flow_fn)
    end
end

############################################################
equations(bgv::BondGraphVertex) = bgv.eqs
equations(::NonParametricJunction) = nothing  # FIXME

numports(bgv::BondGraphVertex) = length(bgv.efforts)
numports(::SourceElement) = 1
numports(::NonParametricJunction) = Inf

############################################################
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

show(io::IO, vertex::T) where {T<:BondGraphVertex} = print(io, "$T{$(numports(vertex))}")
show(io::IO, ::T) where {T<:NonParametricJunction} = print(io, T)

############################################################
struct Port
    name::Symbol
    sys::ODESystem  # port subsys
    parentname::Symbol
    connected::Ref{Bool}
    function Port(name, sys, parentname)
        new(name, sys, parentname, Ref(false))
    end
end
Port(sys::ODESystem, parentname::Symbol) = Port(sys.name, sys, parentname)

name(p::Port) = p.name
system(p::Port; namespaced=true) = namespaced ? ModelingToolkit.renamespace(p.parentname, p.sys) : p.sys
parent(p::Port) = p.parentname

##############################
is_connected(p::Port) = p.connected[]
connect!(p::Port) = p.connected[] = true

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.parentname).$(port.name) $connection_state")
end

############################################################
# Components now define BG elements and junction structures (may rename)
struct Component{T<:BondGraphVertex}
    element::T  # aka component subtype
    name::Symbol
    sys::ODESystem
    ports::Vector{Port}
end
function Component(element::BondGraphVertex; name::Symbol)
    # create MTK System from element definition
    sys = to_system(element; name)
    powerports = getports(sys)
    bg_ports = isempty(powerports) ? Port[] : Port.(powerports, name)
    Component(element, name, sys, bg_ports)
end

##############################

# TODO for sources
""" General port system constructor for Elements --> systems """
function to_system(element::BondElement; name::Symbol)
    # create base system
    sys = System(equations(element), t; name)

    # create ports and add port connections
    port_connection_eqs = Equation[]
    for (i, (e, f)) in enumerate(zip(element.efforts, element.flows))
        # create N port "systems" and extend the user-given MTK System
        port = PowerPort(name=Symbol("port_", i))
        sys = compose(sys, port)

        # add effort/flow connections to newly added port variables (assuming efforts and flows are in the desired order)
        append!(port_connection_eqs, [e ~ port.e, f ~ port.f])
    end
    port_eqs_sys = System(port_connection_eqs, t; name=sys.name)

    return extend(port_eqs_sys, sys)
end

""" Empty system for 0- and 1- Junctions """
function to_system(::NonParametricJunction; name::Symbol)
    System(Equation[], t; name)
end

function getports(sys)
    filter(ModelingToolkit.isconnector, ModelingToolkit.get_systems(sys))
end

############################################################
elementtype(::Component{V}) where {V} = V

show(io::IO, comp::Component{<:BondElement}) = print(io, "$(glyph(comp.element))::$(comp.name)")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.element))")

##############################
system(comp::Component) = comp.sys

# System construction for junction components
# Since the number of ports can change, these systems are composed each time they are called
function system(junc::Component{<:NonParametricJunction})
    # create base system + port subsystems
    basesys = junc.sys
    isempty(junc.ports) && return basesys

    # base system without internal connection equations
    portsys = system.(junc.ports, namespaced=false)
    sys = compose(basesys, portsys)

    # create effort and flow conservation laws
    effort_eqs = junc.element.effort_fn(getefforts(sys))
    flow_eqs = junc.element.flow_fn(getflows(sys))
    extend(System([effort_eqs; flow_eqs], t; name=sys.name), sys)
end

##############################

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
    portsys = PowerPort(; name=Symbol("port_$index"))
    port = Port(portsys, junc.name)
    push!(junc.ports, port) # add new port to junction
    port
end

############################################################################################

struct Bond
    src::Port
    dst::Port
    function Bond(srcport::Port, dstport::Port)
        is_connected(srcport) && error("$srcport already connected")
        is_connected(dstport) && error("$dstport already connected")
        connect!(srcport)
        connect!(dstport)
        new(srcport, dstport)
    end
end
function Bond(srccomp::Component, dstcomp::Component)
    hasfreeport(srccomp) || error("$srccomp has no free ports")
    hasfreeport(dstcomp) || error("$dstcomp has no free ports")
    srcport = nextfreeport(srccomp)
    dstport = nextfreeport(dstcomp)
    Bond(srcport, dstport)
end

ports(b::Bond) = b.src, b.dst
vertices(b::Bond) = parent(b.src), parent(b.dst)

function show(io::IO, b::Bond)
    src, dst = vertices(b)
    print(io, "$src ⇀ $dst")
end

# MTK system connector
function connection_equation(b::Bond)
    srcport_sys = system(b.src)
    dstport_sys = system(b.dst)
    ModelingToolkit.connect(srcport_sys, dstport_sys)
end

############################################################

struct BondGraph
    name::Symbol
    elements::Vector{Component}
    junctions::Vector{Component}
    bonds::Vector{Bond}
    function BondGraph(name, elements=[], junctions=[], bonds=[])
        new(Symbol(name), elements, junctions, bonds)
    end
end
# 3 argument constructor
function BondGraph(name, elements_junctions, bonds)
    elements = filter(x -> x.element isa BondElement, elements_junctions)
    junctions = filter(x -> x.element isa JunctionStructure, elements_junctions)
    BondGraph(Symbol(name), elements, junctions, bonds)
end

# maybe not as the default "show" option (summary print instead)
function show(io::IO, bg::BondGraph)
    print_str = "BondGraph \"$(bg.name)\""
    print_str *= isempty(bg.bonds) ? "" : "\n$(join(bg.bonds,"\n"))"
    print(io, print_str)
end

components(bg::BondGraph) = [bg.elements; bg.junctions]

# constitutive_relations(comp::Component) = equations(system(comp))  # FIXME should be toplevel only


# MTK System converter
function system(bg::BondGraph; simplify=true)
    comps = components(bg)
    subsyss = system.(comps)

    conn_eqns = connection_equation.(bg.bonds)
    basesys = System(conn_eqns, t, name=bg.name)

    sys = compose(basesys, subsyss...)
    simplify ? structural_simplify(sys) : sys
end

############################################################
# TODO CONTINUE FROM HERE
# Graph representation
function graph(bg::BondGraph)
    bg_graph = MetaGraph(
        DiGraph();
        label_type=Symbol,
        vertex_data_type=Component,
        edge_data_type=Bond,
        graph_data=string(bg.name),
        # optional: add weight function and default weight
    )
    for comp in components(bg)
        bg_graph[comp.name] = comp
    end
    for bond in bg.bonds
        src, dst = vertices(bond)
        bg_graph[src.name, dst.name] = bond
    end
    bg_graph
end

############################################################################
using Plots, GraphRecipes
import GraphRecipes: graphplot

function graphplot(bg::BondGraph; kwargs...)
    g = graph(bg)
    graphplot(g;
        title = bg.name,
        names = g.vertex_labels,
        curves = false,
        nodeshape = :rect,
        kwargs...
    )
end
# graphplot(bg)
