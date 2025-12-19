# Base.eltype
eltype(::Type{BondGraph}) = AbstractNode
eltype(::BondGraph) = AbstractNode

# edgetype
edgetype(::Type{BondGraph}) = Graphs.AbstractSimpleEdge{Int}
edgetype(::BondGraph) = Graphs.AbstractSimpleEdge{Int}

# edges
edges(bg::BondGraph) = bg.bonds
ne(bg::BondGraph) = length(bg.bonds)
has_edge(bg::BondGraph, bond::Bond) = any(b -> b === bond, bg.bonds) # strong equality
function has_edge(bg::BondGraph, n1::AbstractNode, n2::AbstractNode)
    has_edge(bg, vertex(n1), vertex(n2))
end
has_edge(bg::BondGraph, s::Int, d::Int) = any(b -> src(b) === s && dst(b) === d, bg.bonds)

# vertices
vertices(bg::BondGraph) = vertex.(bg.nodes)
nv(bg::BondGraph) = length(bg.nodes)
has_vertex(bg::BondGraph, node::AbstractNode) = any(n -> n === node, bg.nodes) # strong equality
has_vertex(bg::BondGraph, v::Int) = 1 <= v <= length(bg.nodes)

# inneighbors, outneighbors
inneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[inneighbors(bg, vertex(n))]
inneighbors(bg::BondGraph, v::Int) = [src(b) for b in bg.bonds if dst(b) == v]
outneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[outneighbors(bg, vertex(n))]
outneighbors(bg::BondGraph, v::Int) = [dst(b) for b in bg.bonds if src(b) == v]
all_neighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[all_neighbors(bg, vertex(n))]

# is_directed
is_directed(::Type{BondGraph}) = true
is_directed(::BondGraph) = true

# zero
zero(::Type{BondGraph}) = BondGraph()
zero(::BondGraph) = BondGraph()

# src, dst
src(b::Bond) = vertex(srcnode(b))
dst(b::Bond) = vertex(dstnode(b))

# Mutations
function add_vertex!(bg::BondGraph, node::AbstractNode)
    has_vertex(bg, node) && return false
    push!(bg.nodes, node)
    set_vertex!(node, nv(bg))
    return true
end

function rem_vertex!(bg::BondGraph, node::AbstractNode)
    has_vertex(bg, node) || return false
    index = vertex(node)
    deleteat!(bg.nodes, index)
    for n in bg.nodes[index:end]
        n.vertex[] -= 1
    end
    return true
end

function add_edge!(bg::BondGraph, srctuple, dsttuple)
    new_bond = Bond(srctuple, dsttuple)
    push!(bg.bonds, new_bond)

    srcnode(new_bond) isa Junction && set_weight!(srctuple..., -1)
    dstnode(new_bond) isa Junction && set_weight!(dsttuple..., +1)

    updateport!(srctuple...)
    updateport!(dsttuple...)

    return new_bond
end

function rem_edge!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    index = findfirst(b -> node1 in b && node2 in b, bg.bonds)
    isnothing(index) && return false # already disconnected

    deleted_bond = bg.bonds[index]
    deleteat!(bg.bonds, index)

    for (node, label) in deleted_bond
        if node isa Junction
            set_weight!(node, label, 0)
        else
            updateport!(node, label)
        end
    end

    return deleted_bond
end
