import Base: show, size
using Graphs, MetaGraphsNext, Symbolics, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

include("ontology.jl")
include("components.jl")
# include("standardlibrary.jl")
# using .Library

############################################################################

struct Bond
    src::PortType
    dst::PortType
    function Bond(src::PortType, dst::PortType)
        src.connected = true
        dst.connected = true
        new(src, dst)
    end
end
function Bond(src_comp::Element, dst_comp::Element)
    hasfreeport(src_comp) || error("$src_comp has no free ports")
    hasfreeport(dst_comp) || error("$dst_comp has no free ports")
    Bond(nextfreeport(src_comp), nextfreeport(dst_comp))
end

# NOTE: I believe domain_connect() is the right way to connect bonds
# see: https://docs.sciml.ai/ModelingToolkit/stable/tutorials/domain_connections/#Special-Connection-Cases-(domain_connect())
conn = domain_connect(port_a, port_b)
conn

############################################################################

# New bond graph structure
mutable struct NewBondGraph{I <: Integer} <: AbstractGraph{I}
    name::AbstractString
    graph::MetaGraph
    sys::ModelingToolkit.AbstractSystem
end
function NewBondGraph(name)
    # creating MetaGraph in a constructor keeps it type stable
    metagraph = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type = Symbol,  # node name
        vertex_data_type = BondGraphVertex,  # node type
        edge_data_type = Bond,  # bond
        graph_data = name  # tag for the whole graph        # TODO add weight function and default weight
    )
    # default "empty" MTK model which is extended with components
    model = ODESystem(Equation[], t; name = Symbol(name))
    new(name, metagraph, model)
end

name(bg::NewBondGraph) = bg.name
size(bg::NewBondGraph) = (nv(bg.graph), ne(bg.graph))
graph(bg::NewBondGraph) = bg.graph

show(io::IO, bg::NewBondGraph) = print(io, "$(name(bg)) BondGraph$(size(bg))")

bg = NewBondGraph("RC Circuit")

############################################################################

function add_node!(bg::NewBondGraph, comp::Element)
    bg.graph[comp.name] = comp
end

function connect!(bg::NewBondGraph, src_comp::Element, dst_comp::Element)
    # TODO assuming single ports, change to allow selecting a specific port
    bg.graph[src_comp.name, dst_comp.name] = Bond(src_comp, dst_comp)
end

############################################################################
############################################################################

@named rc_model = NewBondGraph()

capacitor = Capacitor()
resistor = Resistor()
kvl = ZeroJunction()

add_node!(rc_model, capacitor)
add_node!(rc_model, resistor)
add_node!(rc_model, kvl)

connect!(rc_model, capacitor, kvl)
connect!(rc_model, kvl, resistor)

rc_model

# Graphs.jl
G = graph(rc_model)
incidence_matrix(G)

G[:capacitor]
G.vertex_properties
G.edge_data
G.vertex_labels

############################################################################
using Plots, GraphRecipes
import GraphRecipes: graphplot

function graphplot(bg::NewBondGraph; kwargs...)
    g = graph(bg)
    graphplot(g;
        title = name(bg),
        names = g.vertex_labels,
        curves = false,
        nodeshape = :rect,
        kwargs...
    )
end

graphplot(rc_model)

# TODO CONTINUE FROM HERE
