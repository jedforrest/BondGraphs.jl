import Base: show
using Graphs
using MetaGraphsNext
using Symbolics

# EnergyPair
@variables t
D = Differential(t)

struct EnergyTuple
    p::Num
    q::Num
    e::Num
    f::Num
end
function EnergyTuple(p::Symbol=:p, q::Symbol=:q, e::Symbol=:e, f::Symbol=:f)
    p, q, e, f = @variables $p(t), $q(t), $e(t), $f(t)
    EnergyTuple(p, q, e, f)
end
momentum(et::EnergyTuple) = et.p
position(et::EnergyTuple) = et.q
effort(et::EnergyTuple) = et.e
flow(et::EnergyTuple) = et.f
show(io::IO, et::EnergyTuple) = print(io, "E($(et.p),$(et.q),$(et.e),$(et.f))")

Symbolics.derivative(et::EnergyTuple) = EnergyTuple(et.e, et.f, D(et.e), D(et.f))
power(et::EnergyTuple) = effort(et) * flow(et)
et = EnergyTuple(:l, :x, :u, :v)
power(et)

# Nodes
default_equations(et::EnergyTuple) = [D(et.p) ~ et.e, D(et.q) ~ et.f]
default_equations(et)

abstract type AbstractNode end

struct Component{C} <: AbstractNode
    name::AbstractString
    variables::Vector{EnergyTuple}
    equations::Vector{Equation}
    ports::Vector{Integer}
    function Component(class::Symbol; name::AbstractString="", numports::Integer=1)
        ports = collect(Integer, 1:numports)
        variables = [
            EnergyTuple(Symbol("p$i"), Symbol("q$i"), Symbol("e$i"), Symbol("f$i"))
            for i in 1:numports
        ]
        new{class}(name, variables, Equation[], ports)
    end
end
class(node::AbstractNode) = typeof(node).parameters[1]
numports(node::AbstractNode) = length(node.ports)

show(io::IO, comp::Component) = print(io, "$(class(comp)):{$(join(comp.equations,", "))}")
c = Component(:C)
r = Component(:R; numports=2)
numports(r)

struct Junction{C} <: AbstractNode
    name::AbstractString
    ports::Vector{Integer}
    function Junction(class::Symbol; name::AbstractString="")
        new{class}(name, Integer[])
    end
end
show(io::IO, junc::Junction) = print(io, "$(class(junc))")
Junction(:J0)

# Bonds
struct Bond
    src::Pair{AbstractNode, Integer}  # Component => Port index
    dst::Pair{AbstractNode, Integer}
end
show(io::IO, bond::Bond) = print(io, "<$(src), $(dst)>")

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
bg.graph[:C1] = Component(:C)
bg.graph[:C2] = Component(:C)
bg.graph[:R] = Component(:R, numports=2)
bg.graph[:J0] = Junction(:J0)

# default (linear) equations
default_equations(c::Component) = default_equations.(Ref(c), c.variables)
default_equations(::Component, ::EnergyTuple) = zero(Num)
function default_equations(::Component{:C}, et::EnergyTuple)
    @variables C
    effort(et) ~ C * position(et)
end
function default_equations(::Component{:R}, et::EnergyTuple)
    @variables R
    effort(et) ~ R * flow(et)
end

x = Component(:X)
c = Component(:C)
r = Component(:R, numports=2)
default_equations(x, et)
default_equations(c, et)
default_equations(r)
