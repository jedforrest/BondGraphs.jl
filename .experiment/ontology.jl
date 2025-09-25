# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
import Base: show, size
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: equations

############################################################
# TODO?: https://discourse.julialang.org/t/extracting-kwargs-from-anonymous-function/37350/8

abstract type BondGraphVertex end

abstract type BondElement <: BondGraphVertex end
abstract type StorageElement <: BondElement end
abstract type SourceElement <: BondElement end

"""`R` component"""
struct DissipatorElement <: BondElement
    cr::Any
    parameters::Vector{Num}
    numports::Int
end
function (de::DissipatorElement)(e, f)
    [de.cr(e, f, de.parameters...) ~ 0]
end

"""`C` component"""
struct StaticStorageElement <: StorageElement
    cr::Any
    parameters::Vector{Num}
    numports::Int
end
function (sse::StaticStorageElement)(e, f)
    @variables q(t)
    [sse.cr(e, q, sse.parameters...) ~ 0, D(q) ~ f]
end

"""`I` component"""
struct DynamicStorageElement <: StorageElement
    cr::Vector{Equation}
end

"""`Se` component"""
struct EffortSource <: SourceElement end
"""`Sf` component"""
struct FlowSource <: SourceElement end
"""`SS` component"""
struct SourceSensor <: SourceElement end

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

# Julia Type trees
# using GraphRecipes, Plots
# default(size=(1000, 1000))
# plot(BondGraphVertexClass, method=:tree, fontsize=10, nodeshape=:ellipse)

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

@connector EffortFlow begin
    e(t) = 0., [connect = Effort]
    f(t) = 0., [connect = Flow]
end

struct Port
    parent::Any  # should be Component
    index::Int
    connected::Base.Ref{Bool}
    sys::ODESystem
    function Port(parent::Any, index::Int=1)
        portname = Symbol(parent.name, :_p, index)
        sys = EffortFlow(name=portname)
        new(parent, index, Ref(false), sys)
    end
end
is_connected(p::Port) = p.connected[]
connect!(p::Port) = p.connected[] = true
parent(p::Port) = p.parent
effort(p::Port) = p.sys.e
flow(p::Port) = p.sys.f
vars(p::Port) = effort(p), flow(p)

function show(io::IO, port::Port)
    connection_state = is_connected(port) ? "⬤" : "◯"
    print(io, "⟨$(effort(port)), $(flow(port))⟩ $connection_state")
end

############################################################
# TODO store System in Component definition

# Components now define BG elements and junction structures
struct Component{V<:BondGraphVertex}
    type::V  # rename
    name::Symbol
    ports::Vector{Port}
end
function Component(type::BondElement, name)
    newcomp = Component(type, Symbol(name), Port[])
    for i in 1:numports(type)
        push!(newcomp.ports, Port(newcomp, i))
    end
    newcomp
end
function Component(type::JunctionStructure, name)
    Component(type, Symbol(name), Port[])
end

show(io::IO, comp::Component{<:BondElement}) = print(io, "$(glyph(comp.type))::$(comp.name)")
show(io::IO, comp::Component{<:JunctionStructure}) = print(io, "$(glyph(comp.type))")

subtype(comp::Component) = comp.type
parameters(comp::Component) = comp.type.parameters
variables(comp::Component) = [effort.(comp.ports); flow.(comp.ports)]

hasfreeport(comp::Component) = any(!is_connected, comp.ports)
hasfreeport(::Component{<:NonParametricJunction}) = true

nextfreeport(comp::Component) = first(filter(!is_connected, comp.ports))
function nextfreeport(junc::Component{<:NonParametricJunction})
    # FIXME should only create ports if none are free
    # 0- and 1- junctions have unlimited ports
    # so create a new port if trying to connect
    index = length(junc.ports) + 1
    port = Port(junc, index)
    push!(junc.ports, port)
    port
end

# CRs map (e,f) -> ϕ
function constitutive_relations(comp::Component)
    comptype = subtype(comp)
    es = effort.(comp.ports)
    fs = flow.(comp.ports)
    if numports(comptype) == 1
        return comptype(es[], fs[])
    else
        return comptype(es, fs)
    end
end

# This will eventually become the MTK System converter
function system(elem::Component{<:BondElement})
    eqs = constitutive_relations(elem)
    ODESystem(eqs, t, name=elem.name)
end

# TODO junction should be an MTK connector type
function system(junc::Component{<:JunctionStructure})
    eqs = constitutive_relations(junc)
    ODESystem(eqs, t, name=junc.name)
end

############################################################

struct Bond
    src::Port
    dst::Port
    function Bond(src::Port, dst::Port)
        is_connected(src) && error("$src already connected")
        is_connected(dst) && error("$dst already connected")
        connect!(src)
        connect!(dst)
        new(src, dst)
    end
end
function Bond(src::Component, dst::Component)
    hasfreeport(src) || error("$src has no free ports")
    hasfreeport(dst) || error("$dst has no free ports")
    srcport = nextfreeport(src)
    dstport = nextfreeport(dst)
    Bond(srcport, dstport)
end

vertices(b::Bond) = parent(b.src), parent(b.dst)

function show(io::IO, b::Bond)
    src, dst = vertices(b)
    print(io, "$src ⇀ $dst")
end

# MTK system connector
function connect(b::Bond)
    src_sys = b.src.sys
    dst_sys = b.dst.sys
    ModelingToolkit.connect(src_sys, dst_sys)
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
    elements = filter(x -> x.type isa BondElement, elements_junctions)
    junctions = filter(x -> x.type isa JunctionStructure, elements_junctions)
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


# This will eventually become the MTK System converter
function system(bg::BondGraph; simplify=true)
    comps = components(bg)
    subsyss = system.(comps)

    conn_eqns = connect.(bg.bonds)
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
