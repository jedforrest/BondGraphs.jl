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
        lg.rem_edge!(bg, srcnode(bond), dstnode(bond))
    end
end


function connect!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode)
    srcnode in bg.nodes || error("$srcnode not found in bond graph")
    dstnode in bg.nodes || error("$dstnode not found in bond graph")
    return lg.add_edge!(bg, srcnode, dstnode)
end

function disconnect!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    # rem_edge! removes the bond regardless of the direction of the bond
    deleted_bond = lg.rem_edge!(bg, node1, node2)
    if isnothing(deleted_bond) # if returned nothing, try flipping node1 and node2
        deleted_bond = lg.rem_edge!(bg, node2, node1)
    end
    return deleted_bond
end


function swap!(bg::BondGraph, oldnode::AbstractNode, newnode::AbstractNode)
    numports(oldnode) == numports(newnode) || error("Nodes must have the same number of ports")
    
    srcnodes = lg.inneighbors(bg, oldnode)
    dstnodes = lg.outneighbors(bg, oldnode)

    add_node!(bg, newnode)
    remove_node!(bg, oldnode)

    for src in srcnodes
        connect!(bg, src, newnode)
    end
    for dst in dstnodes
        connect!(bg, newnode, dst)
    end
end


# TODO implement according to https://bondgraphtools.readthedocs.io/en/latest/api.html#BondGraphTools.expose
# function expose!()
    
# end