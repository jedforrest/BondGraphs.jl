function add_nodes!(bg::BondGraph, nodes::AbstractArray{T}) where T <: AbstractNode
    for node in nodes
        add_nodes!(bg, node)
    end
end

function add_nodes!(bg::BondGraph, node::AbstractNode)
    lg.add_vertex!(bg, node) || error("$(typeof(node)) already in model")
end

function connect!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    checkfreeports(bg, node1) || error("$node1 has no free ports")
    checkfreeports(bg, node2) || error("$node2 has no free ports")
    lg.add_edge!(bg, Bond(node1, node2)) || error("Bond already exists between $node1 and $node2")
end