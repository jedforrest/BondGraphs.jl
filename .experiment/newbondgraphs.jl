import Base: show
using Graphs
using MetaGraphsNext
using Symbolics

# EnergyPair
@variables t

struct EnergyPair
    effort::Num
    flow::Num
end
function EnergyPair(effort::Symbol, flow::Symbol)
    e, f = @variables $effort(t), $flow(t)
    EnergyPair(e, f)
end
show(io::IO, ep::EnergyPair) = print(io, "<$(ep.effort), $(ep.flow)>")
power(ep::EnergyPair) = ep.effort * ep.flow
ep = EnergyPair(:u, :v)
power(ep)

# Nodes
abstract type ComponentClass end
struct Capacitance <: ComponentClass end
struct Resistance <: ComponentClass end
struct Inductance <: ComponentClass end
default_equations(ep::EnergyPair, ::ComponentClass ) = power(ep) ~ 0

abstract type AbstractNode end

struct Component <: AbstractNode
    class::ComponentClass
    name::AbstractString
    variables::Vector{EnergyPair}
    equations::Vector{Equation}
    ports::Vector{Integer}
end
function Component(class::ComponentClass; name::AbstractString="", numports::Integer=1)
    ports = collect(Integer, 1:numports)
    variables = [EnergyPair(Symbol("e_$i"), Symbol("f_$i")) for i in 1:numports]
    equations = default_equations.(variables, Ref(class))
    Component(class, name, variables, equations, ports)
end
show(io::IO, comp::Component) = print(io, "$(comp.class)($(comp.equations))")
Component(Capacitance())
Component(Resistance(); numports=2)

struct Bond
    src::Tuple{AbstractNode, Any}
    dst::Tuple{AbstractNode, Any}
end

# New bond graph structure
struct NewBondGraph
    name::AbstractString
    graph::MetaGraph
    function NewBondGraph(name::AbstractString, graph::Graphs.AbstractGraph)
        # creating MetaGraph in constructor keeps it type stable
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
bg.graph[:C1] = Component(Capacitance())
bg.graph[:C2] = Component(Capacitance())
bg.graph[:R] = Component(Capacitance(), numports=2)
