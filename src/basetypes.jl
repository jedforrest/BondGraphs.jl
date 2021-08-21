abstract type AbstractNode end

struct Component <: AbstractNode
    metamodel::Symbol
    name::AbstractString
    maxports::Int
end
Component(metamodel::Symbol, name::String=string(metamodel)) = Component(metamodel, name, 1)

struct Junction <: AbstractNode
    metamodel::Symbol
end

struct Bond <: lg.AbstractSimpleEdge{Integer}
    src::AbstractNode
    dst::AbstractNode
end

struct BondGraph <: lg.AbstractGraph{Int64}
    metamodel::Symbol
    name::AbstractString
    nodes::Vector{T} where T <:AbstractNode
    bonds::Vector{Bond}
end
BondGraph(metamodel::Symbol=:BG; name::String="bg") = BondGraph(metamodel, name, AbstractNode[], Bond[])

# Indexing
find_index(bg::BondGraph, node::AbstractNode) = findfirst(x -> x == node, bg.nodes)
find_index(bg::BondGraph, bond::Bond) = findfirst(x -> x == bond, bg.bonds)

# Ports
maxports(n::Component) = n.maxports
maxports(n::Junction) = Inf
checkfreeports(bg::BondGraph, n::AbstractNode) = length(lg.all_neighbors(bg, find_index(bg, n))) < maxports(n)

# I/O
show(io::IO, node::Component) = print(io, "$(node.metamodel):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.metamodel)")
show(io::IO, b::Bond) = print(io, "Bond $(b.src) ⇀ $(b.dst)")
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.metamodel):$(bg.name) ($(lg.nv(bg)) Nodes, $(lg.ne(bg)) Bonds)")
