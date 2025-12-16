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
    e(t) = 0., [connect = Effort]
    f(t) = 0., [connect = Flow]
end

############################################################
getefforts(sys::ODESystem) = [e for e in unknowns(sys) if getconnect(e) == Effort]
getflows(sys::ODESystem) = [f for f in unknowns(sys) if getconnect(f) == Flow]

""" General port constructor for Bond Elements """
function add_port_connections(sys::ODESystem, numports::Int=1)
    # get effort and flow vars from user-given System
    efforts = getefforts(sys)
    flows = getflows(sys)

    # check validity
    (length(efforts) == length(flows) == numports) || error("Number of efforts, flows, and ports don't match ($numports)")

    port_connection_eqs = Equation[]
    for i in 1:numports
        # create N port "systems" and extend the user-given MTK System
        powerport = PowerPort(name=Symbol("port_", i))
        sys = compose(sys, powerport)

        # add effort/flow connections to newly added port variables (assuming efforts and flows are in the desired order)
        append!(port_connection_eqs, [efforts[i] ~ powerport.e, flows[i] ~ powerport.f])
    end
    port_eqs_sys = System(port_connection_eqs, t; name=sys.name)

    return extend(port_eqs_sys, sys)
end

############################################################
# TODO instead of storing a model or sys, this can store the symbolic equations?
# and label the efforts and flows (and state vars) explicitly

# NOTE: SDESystems do not yet work with @mtkmodel - therefore must use Non-DSL approach
"""`R` component"""
struct DissipatorElement <: BondElement
    sys::ODESystem
    numports::Int
    function DissipatorElement(sys::ODESystem, numports::Int=1)
        sys = add_port_connections(sys, numports)
        new(sys, numports)
    end
end

"""`C` component"""
struct StaticStorageElement <: StorageElement
    sys::ODESystem
    numports::Int
    function StaticStorageElement(sys::ODESystem, numports::Int=1)
        sys = add_port_connections(sys, numports)
        new(sys, numports)
    end
end

"""`I` component"""
struct DynamicStorageElement <: StorageElement
    sys::ODESystem
    numports::Int
    function DynamicStorageElement(sys::ODESystem, numports::Int=1)
        sys = add_port_connections(sys, numports)
        new(sys, numports)
    end
end

# TODO make defaults
res
eqs = equations(res)
sys3 = System(equations(res), t; name=res.name)

V, I = unknowns(res)
R, = ModelingToolkit.parameters(res)

res == sys3

rcomp.element
@variables y(t)

@named r = Library.Resistor()
Library.Resistor

############################################################
"""`Se` component"""
struct EffortSource <: SourceElement
    sys::ODESystem
    function StaticStorageElement(sys::ODESystem)
        sys = add_port_connections(sys, 1)
        new(sys)
    end
end

"""`Sf` component"""
struct FlowSource <: SourceElement
    sys::ODESystem
    function FlowSource(sys::ODESystem)
        sys = add_port_connections(sys, 1)
        new(sys)
    end
end

"""`SS` component"""
struct SourceSensor <: SourceElement
    sys::ODESystem
    function SourceSensor(sys::ODESystem)
        sys = add_port_connections(sys, 1)
        new(sys)
    end
end

############################################################
abstract type JunctionStructure <: BondGraphVertex end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

"""`TF` component"""
struct Transformer <: ParametricJunction
    sys::ODESystem
    numports::Int
    function Transformer(sys::ODESystem, numports::Int=2)
        numports >= 2 || error("Transformer must have at least 2 ports ($N ports given)")
        sys = add_port_connections(sys, numports)
        new(sys, numports)
    end
end
"""`GY` component"""
struct Gyrator <: ParametricJunction
    sys::ODESystem
    numports::Int
    function Gyrator(sys::ODESystem, numports::Int=2)
        numports >= 2 || error("Gyrator must have at least 2 ports ($N ports given)")
        sys = add_port_connections(sys, numports)
        new(sys, numports)
    end
end

############################################################

# 0- and 1- junctions start with "empty" systems
"""`0`-junction"""
struct EqualEffort <: NonParametricJunction
    sys::ODESystem
    function EqualEffort(; name::Symbol=:j0)
        new(System(Equation[], t; name=name))
    end
end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction
    sys::ODESystem
    function EqualFlow(; name::Symbol=:j1)
        new(System(Equation[], t; name=name))
    end
end

# easy way to construct equtions for given effort/flow inputs
function (::EqualEffort)(e, f)
    e_eqs = [e[1] ~ ei for ei in e[2:end]]
    f_eqs = sum(f) ~ 0
    [e_eqs; f_eqs]
end
function (::EqualFlow)(e, f)
    f_eqs = [f[1] ~ fi for fi in f[2:end]]
    e_eqs = sum(e) ~ 0
    [f_eqs; e_eqs]
end

############################################################
equations(bgv::BondGraphVertex) = equations(bgv.sys)

numports(bgv::BondGraphVertex) = bgv.numports
numports(::SourceElement) = 1
numports(::NonParametricJunction) = Inf

function getports(be::BondElement)
    filter(ModelingToolkit.isconnector, ModelingToolkit.get_systems(be.sys))
end

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
struct Component{V<:BondGraphVertex}
    element::V  # aka component subtype
    name::Symbol
    ports::Vector{Port}
end
function Component(element::BondElement; name::Symbol)
    bg_ports = Port.(getports(element), name)
    Component(element, name, bg_ports)
end
function Component(element::JunctionStructure; name::Symbol)
    Component(element, name, Port[])
end

############################################################
elementtype(::Component{V}) where {V} = V

show(io::IO, comp::Component{<:BondElement}) = print(io, "$(glyph(comp.element))::$(comp.name)")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.element))")

system(comp::Component) = comp.element.sys

##############################
"""System construction for junction components. These systems are composed when called."""
function system(junc::Component{<:NonParametricJunction})
    # create base system + port subsystems
    basesys = junc.element.sys
    if isempty(junc.ports)
        return basesys
    end
    # base system without internal connection equations
    portsys = system.(junc.ports, namespaced=false)
    sys = compose(basesys, portsys)

    # TODO cleanup this code (its not very clear)
    e, f = getefforts(sys), getflows(sys)
    inner_connection_eqs = junc.element(e, f)
    extend(System(inner_connection_eqs, t; name=sys.name), sys)
end

##############################

parameters(comp::Component) = ModelingToolkit.parameters(system(comp))
variables(comp::Component) = ModelingToolkit.get_unknowns(system(comp))  # FIXME should be toplevel only

# TODO dispatch based on element subtype
efforts(comp::Component) = getefforts(system(comp))
flows(comp::Component) = getflows(system(comp))

constitutive_relations(comp::Component) = equations(system(comp))  # FIXME should be toplevel only

############################################################
hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, comp.ports))
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    portsys = PowerPort(; name=Symbol("port_$index"))
    port = Port(portsys, junc.name)
    push!(junc.ports, port)
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

# MTK System converter
function system(bg::BondGraph; simplify=true)
    comps = components(bg)
    subsyss = system.(comps)

    conn_eqns = connection_equation.(bg.bonds)
    basesys = ODESystem(conn_eqns, t, name=bg.name)

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
