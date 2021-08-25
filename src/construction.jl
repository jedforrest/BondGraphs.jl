function add_node!(bg::BondGraph, nodes::AbstractArray{T}) where T <: AbstractNode
    for node in nodes
        add_node!(bg, node)
    end
end

function add_node!(bg::BondGraph, node::AbstractNode)
    lg.add_vertex!(bg, node) || error("$(typeof(node)) already in model")
end


function remove_node!(bg::BondGraph, nodes::AbstractArray{T}) where T <: AbstractNode
    for node in nodes
        remove_node!(bg, node)
    end
end

function remove_node!(bg::BondGraph, node::AbstractNode)
    lg.rem_vertex!(bg, node) || error("$(typeof(node)) not in model")
    for bond in filter(bond -> node in bond, bg.bonds)
        lg.rem_edge!(bg, bond)
    end
end


function connect!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    srcnode in bg.nodes || error("$srcnode not found in bong graph")
    dstnode in bg.nodes || error("$dstnode not found in bong graph")
    checkfreeports(bg, srcnode) || error("$srcnode has no free ports")
    checkfreeports(bg, dstnode) || error("$dstnode has no free ports")
    lg.add_edge!(bg, Bond(srcnode, dstnode)) || error("Bond already exists between $node1 and $node2")
end

function disconnect!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    # Try removing an edge in both directions
    lg.rem_edge!(bg, Bond(node1, node2)) && return
    lg.rem_edge!(bg, Bond(node2, node1)) && return
end


function swap!(bg::BondGraph, oldnode::AbstractNode, newnode::AbstractNode)
    maxports(oldnode) == maxports(newnode) || error("Nodes must have the same number of ports")
    
    srcnodes = lg.inneighbors(bg, oldnode)
    dstnodes = lg.outneighbors(bg, oldnode)

    add_node!(bg, newnode)
    remove_node!(bg, oldnode)

    for dst in dstnodes
        connect!(bg, newnode, dst)
    end
    for src in srcnodes
        connect!(bg, src, newnode)
    end
end


# TODO implement according to https://bondgraphtools.readthedocs.io/en/latest/api.html#BondGraphTools.expose
# function expose!()
    
# end