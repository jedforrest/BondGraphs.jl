# Base.eltype
eltype(::Type{BondGraph}) = AbstractNode
eltype(bg::BondGraph) = AbstractNode

# edgetype
lg.edgetype(::Type{BondGraph}) = lg.AbstractSimpleEdge{Integer}
lg.edgetype(bg::BondGraph) = lg.AbstractSimpleEdge{Integer}

# edges
lg.edges(bg::BondGraph) = bg.bonds
lg.ne(bg::BondGraph) = length(bg.bonds)
lg.has_edge(bg::BondGraph, b::Bond) = b in bg.bonds
lg.has_edge(bg::BondGraph, n1::AbstractNode, n2::AbstractNode) = lg.has_edge(bg, vertex(n1), vertex(n2))
lg.has_edge(bg::BondGraph, s::Int, d::Int) = 
    any(b -> lg.src(b) == s && lg.dst(b) == d, bg.bonds)

# vertices
lg.vertices(bg::BondGraph) = vertex.(bg.nodes)
lg.nv(bg::BondGraph) = length(bg.nodes)
lg.has_vertex(bg::BondGraph, v::Int) = v <= length(bg.nodes)
lg.has_vertex(bg::BondGraph, node::AbstractNode) = node in bg.nodes

# inneighbors, outneighbors
lg.inneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[lg.inneighbors(bg, vertex(n))]
lg.inneighbors(bg::BondGraph, v::Int) = [lg.src(b) for b in bg.bonds if lg.dst(b) == v]
lg.outneighbors(bg::BondGraph, n::AbstractNode) = bg.nodes[lg.outneighbors(bg, vertex(n))]
lg.outneighbors(bg::BondGraph, v::Int) = [lg.dst(b) for b in bg.bonds if lg.src(b) == v]

# is_directed
lg.is_directed(::Type{BondGraph}) = true
lg.is_directed(bg::BondGraph) = true

# zero
lg.zero(::Type{BondGraph}) = BondGraph()
lg.zero(bg::BondGraph) = BondGraph()

# src, dst
lg.src(b::Bond) = vertex(srcnode(b))
lg.dst(b::Bond) = vertex(dstnode(b))

# weights
# TODO

# Mutations
function lg.add_vertex!(bg::BondGraph, node::AbstractNode)
    lg.has_vertex(bg, node) && return false
    push!(bg.nodes, node)
    set_vertex!(node, lg.nv(bg))
    return true
end

function lg.rem_vertex!(bg::BondGraph, node::AbstractNode)
    lg.has_vertex(bg, node) || return false
    index = vertex(node)
    deleteat!(bg.nodes, index)
    for n in bg.nodes[index:end]
        n.vertex[] -= 1
    end
    return true
end

function lg.add_edge!(bg::BondGraph, srcport::Port, dstport::Port)
    srcnode = srcport.node
    dstnode = dstport.node
    lg.has_edge(bg, srcnode, dstnode) && error("Bond already exists between $srcnode and $dstnode")
    lg.has_edge(bg, dstnode, srcnode) && error("Bond already exists between $dstnode and $srcnode")
    new_bond = Bond(srcport, dstport)
    push!(bg.bonds, new_bond)
    updateport!(srcnode, srcport.index)
    updateport!(dstnode, dstport.index)
    return new_bond
end
function lg.add_edge!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    lg.add_edge!(bg, Port(srcnode), Port(dstnode))
end

function lg.rem_edge!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    lg.has_edge(bg, srcnode, dstnode) || lg.has_edge(bg, dstnode, srcnode) || return
    index = findfirst(b -> (srcnode in b && dstnode in b), bg.bonds)
    deleted_bond = bg.bonds[index]
    deleteat!(bg.bonds, index)
    updateport!(srcnode, deleted_bond.src.index)
    updateport!(dstnode, deleted_bond.dst.index)
    return deleted_bond
end
