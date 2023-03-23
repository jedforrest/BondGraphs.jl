# Base.eltype
eltype(::Type{BondGraph}) = AbstractNode
eltype(::BondGraph) = AbstractNode

# edgetype
g.edgetype(::Type{BondGraph}) = g.AbstractSimpleEdge{Integer}
g.edgetype(::BondGraph) = g.AbstractSimpleEdge{Integer}

# edges
g.edges(bg::BondGraph) = bg.bonds
g.ne(bg::BondGraph) = length(bg.bonds)
g.has_edge(bg::BondGraph, bond::Bond) =  any(b -> b === bond, bg.bonds) # strong equality
g.has_edge(bg::BondGraph, n1::AbstractNode, n2::AbstractNode) = g.has_edge(bg, vertex(n1), vertex(n2))
g.has_edge(bg::BondGraph, s::Int, d::Int) =
    any(b -> g.src(b) === s && g.dst(b) === d, bg.bonds)

# vertices
g.vertices(bg::BondGraph) = vertex.(bg.nodes)
g.nv(bg::BondGraph) = length(bg.nodes)
g.has_vertex(bg::BondGraph, node::AbstractNode) = any(n -> n === node, bg.nodes) # strong equality
g.has_vertex(bg::BondGraph, v::Int) = 1 <= v <= length(bg.nodes)

# inneighbors, outneighbors
g.inneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.inneighbors(bg, vertex(n))]
g.inneighbors(bg::BondGraph, v::Int) = [g.src(b) for b in bg.bonds if g.dst(b) == v]
g.outneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.outneighbors(bg, vertex(n))]
g.outneighbors(bg::BondGraph, v::Int) = [g.dst(b) for b in bg.bonds if g.src(b) == v]
g.all_neighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.all_neighbors(bg, vertex(n))]

# is_directed
g.is_directed(::Type{BondGraph}) = true
g.is_directed(::BondGraph) = true

# zero
g.zero(::Type{BondGraph}) = BondGraph()
g.zero(::BondGraph) = BondGraph()

# src, dst
g.src(b::Bond) = vertex(srcnode(b))
g.dst(b::Bond) = vertex(dstnode(b))

# weights
# TODO

# Mutations
function g.add_vertex!(bg::BondGraph, node::AbstractNode)
    g.has_vertex(bg, node) && return false
    push!(bg.nodes, node)
    set_vertex!(node, g.nv(bg))
    return true
end

function g.rem_vertex!(bg::BondGraph, node::AbstractNode)
    g.has_vertex(bg, node) || return false
    index = vertex(node)
    deleteat!(bg.nodes, index)
    for n in bg.nodes[index:end]
        n.vertex[] -= 1
    end
    return true
end

function g.add_edge!(bg::BondGraph, srctuple, dsttuple)
    new_bond = Bond(srctuple, dsttuple)
    push!(bg.bonds, new_bond)

    srcnode(new_bond) isa Junction && set_weight!(srctuple..., -1)
    dstnode(new_bond) isa Junction && set_weight!(dsttuple..., +1)

    updateport!(srctuple...)
    updateport!(dsttuple...)

    return new_bond
end

function g.rem_edge!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
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
