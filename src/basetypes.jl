abstract type AbstractNode end

struct Component <: AbstractNode
    metamodel::Symbol
    name::AbstractString
    maxports::Int
    vertex::RefValue{Int}
    function Component(m::Symbol, n::AbstractString, mp::Int, v::Int)
        new(m, n, mp, Ref(v))
    end
end
Component(metamodel::Symbol; name::String=string(metamodel), maxports::Int=1, vertex::Int=0) = 
    Component(metamodel, name, maxports, vertex)

struct Junction <: AbstractNode
    metamodel::Symbol
    vertex::RefValue{Int}
    Junction(m::Symbol; v::Int=0) = new(m, Ref(v))
end

struct Bond <: lg.AbstractSimpleEdge{Int}
    srcnode::AbstractNode
    dstnode::AbstractNode
end

struct BondGraph <: lg.AbstractGraph{Int64}
    metamodel::Symbol
    name::AbstractString
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(metamodel::Symbol=:BG; name::String="bg") = BondGraph(metamodel, name, AbstractNode[], Bond[])

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Ports
maxports(n::Component) = n.maxports
maxports(n::Junction) = Inf
checkfreeports(bg::BondGraph, n::AbstractNode) = length(lg.all_neighbors(bg, vertex(n))) < maxports(n)

# Comparisons
in(n::AbstractNode, b::Bond) = n == b.srcnode || n == b.dstnode

# I/O
show(io::IO, node::Component) = print(io, "$(node.metamodel):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.metamodel)")
show(io::IO, b::Bond) = print(io, "Bond $(b.srcnode) ⇀ $(b.dstnode)")
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.metamodel):$(bg.name) ($(lg.nv(bg)) Nodes, $(lg.ne(bg)) Bonds)")
