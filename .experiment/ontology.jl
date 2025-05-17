# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4
using Graphs, Symbolics

include("energypair.jl")

abstract type BondGraphElement end
abstract type StorageElement <: BondGraphElement end
abstract type DissipatorElement <: BondGraphElement end
abstract type SourceElement <: BondGraphElement end

abstract type JunctionStructure end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end


struct Port
    name::Symbol
    variables::EnergyPair
    parent::Symbol
    connected::Bool
end
Port(name, parent) = Port(Symbol(name), EnergyPair(), Symbol(parent), false)

struct Bond
    src::Port
    dst::Port
end


abstract type AbstractNode end

struct Component{T<:BondGraphElement} <: AbstractNode
    name::Symbol
    ports::Vector{Port}
    equations::Vector{Equation}
    parameters::Vector{Num}
    function Component{T<:BondGraphElement}(
            numports::Integer=1;
            name,
            equations=Equation[],
            parameters=Num[]
        )
        ports = map(i -> Port("p$i", name), 1:numports)
        new{T}(name, ports, equations, parameters)
    end
end

struct Junction{T<:JunctionStructure} <: AbstractNode
    name::Symbol
    ports::Vector{Port}
end


# New bond graph structure
struct NewBondGraph
    name::Symbol
    graph::MetaGraph
    function NewBondGraph(graph::AbstractGraph; name::Symbol)
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
NewBondGraph(; name::Symbol=:BG) = NewBondGraph(Graph(); name)
