# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
import Base: show, size
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: equations, Model, System
using DifferentialEquations

include("./standardlibrary.jl")
using .Library

############################################################
abstract type AbstractPowerVariableType end

struct EffortVar <: AbstractPowerVariableType end
struct FlowVar <: AbstractPowerVariableType end
struct MomentumVar <: AbstractPowerVariableType end
struct DisplacementVar <: AbstractPowerVariableType end

############################################################
abstract type BondGraphVertex end

abstract type BondElement <: BondGraphVertex end
abstract type StorageElement <: BondElement end
abstract type SourceElement <: BondElement end

# NOTE: SDESystems do not yet work with @mtkmodel - therefore must use Non-DSL approach

const Effort = ModelingToolkit.Equality

@connector PowerPort begin
    e(t) = 0., [connect = Effort]
    f(t) = 0., [connect = Flow]
end

# ############################################################
# # TODO move to Library
# @component function Res(; name)
#     vars = @variables begin
#         V(t), [connect = Effort]
#         I(t), [connect = Flow]
#     end
#     ps = @parameters begin
#         R = 2
#     end
#     eqs = [V ~ R * I]
#     return System(eqs, t, vars, ps; name)
# end
# ############################################################


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
@named res = Library.Resistor()
DissipatorElement(res)

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
# TODO CONTINUE FROM HERE
abstract type JunctionStructure <: BondGraphVertex end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

"""`TF` component"""
struct Transformer <: ParametricJunction
    cr::Any
    numports::Int
end
"""`GY` component"""
struct Gyrator <: ParametricJunction
    cr::Any
    numports::Int
end

"""`0`-junction"""
struct EqualEffort <: NonParametricJunction end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction end

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

numports(v::BondGraphVertex) = v.numports
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
const Effort = ModelingToolkit.Equality

@connector PortVars begin
    e(t), [connect = Effort]
    f(t), [connect = Flow]
end

struct Port
    name::Symbol
    parent::Any  # should be Component
    connected::Ref{Bool}
    function Port(name, parent)
        new(Symbol(name), parent, Ref(false))
    end
end
is_connected(p::Port) = p.connected[]
connect!(p::Port) = p.connected[] = true
parent(p::Port) = p.parent
system(p::Port) = getproperty(system(p.parent), p.name)


function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "$(port.parent.name).$(port.name) $connection_state")
end

############################################################
# TODO? Numports as a parametric type?

# Components now define BG elements and junction structures
struct Component{V<:BondGraphVertex}
    subtype::V  # rename
    name::Symbol
    ports::Vector{Port}
end
function Component(type::BondElement; name)
    newcomp = Component(type, Symbol(name), Port[])
    for i in 1:numports(type)
        push!(newcomp.ports, Port("_$i", newcomp))
    end
    newcomp
end
function Component(comptype::JunctionStructure; name::Symbol)
    Component(comptype, name, Port[])
end

show(io::IO, comp::Component{<:BondElement}) = print(io, "$(glyph(comp.subtype))::$(comp.name)")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.subtype))")

subtype(comp::Component) = comp.subtype
parameters(comp::Component) = comp.subtype.parameters
variables(comp::Component) = [effort.(comp.ports); flow.(comp.ports)]

hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, comp.ports))
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    port = Port("_$index", junc)
    push!(junc.ports, port)
    port
end

############################################################


# CRs map (e,f) -> ϕ
# TODO is this function useful?
function constitutive_relations(comp::Component)
    equations(toggle_namespacing(system(comp), false))  # FIXME
end

# This will eventually become the MTK System converter
function system(comp::Component)
    comptype = subtype(comp)

    port_sys = [PortVars(name=port.name) for port in comp.ports]
    es = [ps.e for ps in port_sys]
    fs = [ps.f for ps in port_sys]
    if numports(comptype) == 1
        eqs = comptype(es[], fs[])
    else
        eqs = comptype(es, fs)
    end
    compose(ODESystem(eqs, t; name=comp.name), port_sys)
end

############################################################

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

# MTK system connector (rename)
function connect_equation(b::Bond)
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
    elements = filter(x -> x.subtype isa BondElement, elements_junctions)
    junctions = filter(x -> x.subtype isa JunctionStructure, elements_junctions)
    BondGraph(Symbol(name), elements, junctions, bonds)
end
# 2 argument constructor
function BondGraph(name, bonds::Vector{Bond})
    src_dsts = reduce(vcat, collect.(vertices.(bonds)))
    BondGraph(Symbol(name), unique(src_dsts), bonds)
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

    conn_eqns = connect_equation.(bg.bonds)
    basesys = ODESystem(conn_eqns, t, name=bg.name)

    sys = compose(basesys, subsyss...)
    simplify ? structural_simplify(sys) : sys
end


############################################################
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
