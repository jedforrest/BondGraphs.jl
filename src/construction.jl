function add_node!(bg::BondGraph, nodes)
    for node in nodes
        add_node!(bg, node)
    end
end

function add_node!(bg::BondGraph, node::AbstractNode)
    g.add_vertex!(bg, node) || @warn "Node '$(name(node))' already in model"
end


function remove_node!(bg::BondGraph, nodes)
    for node in nodes
        remove_node!(bg, node)
    end
end

function remove_node!(bg::BondGraph, node::AbstractNode)
    g.rem_vertex!(bg, node) || @warn "Node '$(name(node))' not in model"
    for bond in filter(bond -> node in bond, bg.bonds)
        g.rem_edge!(bg, srcnode(bond), dstnode(bond))
    end
end


function connect!(bg::BondGraph, srcnode::AbstractNode, dstnode::AbstractNode;
        srcportindex=nextfreeport(srcnode), dstportindex=nextfreeport(dstnode))
    srcnode in bg.nodes || error("$srcnode not found in bond graph")
    dstnode in bg.nodes || error("$dstnode not found in bond graph")
    srcport = Port(srcnode, srcportindex)
    dstport = Port(dstnode, dstportindex)
    return g.add_edge!(bg, srcport, dstport)
end

function disconnect!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    # rem_edge! removes the bond regardless of the direction of the bond
    return g.rem_edge!(bg, node1, node2)
end

# TODO
# Flip bond function

function swap!(bg::BondGraph, oldnode::AbstractNode, newnode::AbstractNode)
    _check_port_number(oldnode,newnode)

    # may be a redundant check
    if !g.has_vertex(bg, newnode)
        add_node!(bg, newnode)
    end

    srcnodes = g.inneighbors(bg, oldnode)
    dstnodes = g.outneighbors(bg, oldnode)
    remove_node!(bg, oldnode)

    for src in srcnodes
        connect!(bg, src, newnode)
    end
    for dst in dstnodes
        connect!(bg, newnode, dst)
    end
end

_check_port_number(oldnode::AbstractNode, newnode::AbstractNode) =
    numports(newnode) >= numports(oldnode) || error("New node must have a greater or equal number of ports to the old node")
_check_port_number(oldnode::AbstractNode, newnode::Junction) = true


# TODO implement according to https://bondgraphtools.readthedocs.io/en/latest/api.html#BondGraphTools.expose
# function expose!()

# end


# Inserts an AbstractNode between two connected (bonded) nodes
# The direction of the original bond is preserved by this action
function insert_node!(bg::BondGraph, bond::Bond, newnode::AbstractNode)
    src = srcnode(bond)
    dst = dstnode(bond)

    disconnect!(bg, src, dst)

    try
        add_node!(bg, newnode)
        connect!(bg, src, newnode)
        connect!(bg, newnode, dst)
    catch e
        # if connection fails, reconnect original bond
        disconnect!(bg, src, newnode)
        disconnect!(bg, newnode, dst)
        connect!(bg, src, dst)
        error(e)
    end
end
function insert_node!(bg::BondGraph, tuple::Tuple, newnode::AbstractNode)
    bonds = getbonds(bg, tuple)
    isempty(bonds) && error("$(tuple[1]) and $(tuple[2]) are not connected")
    insert_node!(bg, bonds[1], newnode)
end


function merge_nodes!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode; junction=EqualEffort())
    node1.type == node2.type || error("$(node1.name) must be the same type as $(node2.name)")

    # node1 taken as the node to keep
    for nb in g.all_neighbors(bg, node1)
        junc = deepcopy(junction)
        bond = getbonds(bg, node1, nb)[1]
        insert_node!(bg, bond, junc)
        swap!(bg, node2, junc)
    end
end
function merge_nodes!(bg::BondGraph, node1::Junction, node2::Junction)
    # node1 taken as the node to keep
    # remove conflicting connections between junctions if they exist
    disconnect!(bg, node1, node2)
    shared_neighbors = intersect(g.all_neighbors(bg, node1), g.all_neighbors(bg, node2))
    for shared_neighbor in shared_neighbors
        disconnect!(bg, node2, shared_neighbor)
    end
    swap!(bg, node2, node1)
end


function simplify_junctions!(bg::BondGraph; remove_redundant=true, squash_identical=true)
    junctions = filter(n -> n isa Junction, bg.nodes)

    # Removes junctions with 2 or less connected ports
    if remove_redundant
        for j in junctions
            n_nbrs = length(g.all_neighbors(bg, j))
            if n_nbrs == 2
                #srcnode = g.inneighbors(bg, j)[1]
                #dstnode = g.outneighbors(bg, j)[1]
                node1, node2 = g.all_neighbors(bg, j)
                remove_node!(bg, j)
                # bond direction may not be preserved here
                connect!(bg, node1, node2)
            elseif n_nbrs < 2
                remove_node!(bg, j)
            end
        end
    end

    # Squashes identical copies of the same junction type into one junction
    if squash_identical
        for j in junctions, nbr in g.all_neighbors(bg, j)
            g.has_vertex(bg, j) || continue # in case j was removed
            if type(j) == type(nbr)
                merge_nodes!(bg, j, nbr)
            end
        end
    end
    bg
end
