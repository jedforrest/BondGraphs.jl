# Base.eltype
eltype(::Type{BondGraph}) = AbstractNode
eltype(bg::BondGraph) = AbstractNode

# edgetype
g.edgetype(::Type{BondGraph}) = g.AbstractSimpleEdge{Integer}
g.edgetype(bg::BondGraph) = g.AbstractSimpleEdge{Integer}

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
g.has_vertex(bg::BondGraph, v::Int) = v <= length(bg.nodes)

# inneighbors, outneighbors
g.inneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.inneighbors(bg, vertex(n))]
g.inneighbors(bg::BondGraph, v::Int) = [g.src(b) for b in bg.bonds if g.dst(b) == v]
g.outneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.outneighbors(bg, vertex(n))]
g.outneighbors(bg::BondGraph, v::Int) = [g.dst(b) for b in bg.bonds if g.src(b) == v]
g.all_neighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[g.all_neighbors(bg, vertex(n))]

# is_directed
g.is_directed(::Type{BondGraph}) = true
g.is_directed(bg::BondGraph) = true

# zero
g.zero(::Type{BondGraph}) = BondGraph()
g.zero(bg::BondGraph) = BondGraph()

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

function g.add_edge!(bg::BondGraph, srcport::Port, dstport::Port)
    srcnode = srcport.node
    dstnode = dstport.node
    g.has_edge(bg, srcnode, dstnode) && error("Bond already exists between $srcnode and $dstnode")
    g.has_edge(bg, dstnode, srcnode) && error("Bond already exists between $dstnode and $srcnode")
    new_bond = Bond(srcport, dstport)
    push!(bg.bonds, new_bond)
    updateport!(srcnode, srcport.index, Out)
    updateport!(dstnode, dstport.index, In)
    return new_bond
end
function g.add_edge!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    g.add_edge!(bg, Port(srcnode), Port(dstnode))
end

function g.rem_edge!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    g.has_edge(bg, srcnode, dstnode) || g.has_edge(bg, dstnode, srcnode) || return
    index = findfirst(b -> srcnode in b && dstnode in b, bg.bonds)
    deleted_bond = bg.bonds[index]
    deleteat!(bg.bonds, index)
    updateport!(srcnode, deleted_bond.src.index, Free)
    updateport!(dstnode, deleted_bond.dst.index, Free)
    return deleted_bond
end
