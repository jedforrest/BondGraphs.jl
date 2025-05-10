import Base: show
using Graphs
using MetaGraphsNext

include("newcomponents.jl")

# Bonds
struct Bond
    src::Pair{AbstractNode, Integer}  # Component => Port index
    dst::Pair{AbstractNode, Integer}
end

# New bond graph structure
struct NewBondGraph
    name::AbstractString
    graph::MetaGraph
    function NewBondGraph(name::AbstractString, graph::Graphs.AbstractGraph)
        # creating MetaGraph in a constructor keeps it type stable
        metagraph = MetaGraph(
            graph;  # underlying graph structure
            label_type=Symbol,  # node name
            vertex_data_type=AbstractNode,  # node type
            edge_data_type=Bond,  # bond
            graph_data=name,  # tag for the whole graph
        )
        return new(name, metagraph)
    end
end
NewBondGraph(name::AbstractString="New BG") = NewBondGraph(name, Graph())
show(io::IO, bg::NewBondGraph) = print(io, "BondGraph($(bg.name))")

bg = NewBondGraph()
bg.graph[:C1] = Component(:C)
bg.graph[:C2] = Component(:C)
bg.graph[:R] = Component(:R, 2)
bg.graph[:J0] = Junction(:J0)

x = Component(:X)
c = Component(:C)
r = Component(:R, 2)
