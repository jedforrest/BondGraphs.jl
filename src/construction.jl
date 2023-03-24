"""
    add_node!(bg::BondGraph, nodes)

Add a node to a bond graph `bg`. Can add a single node or list of nodes.
"""
function add_node!(bg::BondGraph, nodes)
    for node in nodes
        add_node!(bg, node)
    end
end

function add_node!(bg::BondGraph, node::AbstractNode)
    g.add_vertex!(bg, node) || @warn "Node '$(name(node))' already in model"
end

"""
    remove_node!(bg::BondGraph, nodes)

Remove a node to a bond graph `bg`. Can remove a single node or list of nodes.
"""
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

"""
    connect!(bg::BondGraph, source_node, destination_node)
    connect!(bg::BondGraph, (source_node, port_label), (destination_node, port_label))

Connect two components together in the same bond graph. The bond direction is always from
`source_node` to `destination_node`. The port index of `source_node` and `destination_node`
can be optionally set.
"""
function connect!(bg::BondGraph, src, dst)
    (srcnode, srcport) = port_info(src)
    (dstnode, dstport) = port_info(dst)

    srcnode in nodes(bg) || error("$srcnode not found in bond graph")
    dstnode in nodes(bg) || error("$dstnode not found in bond graph")
    isnothing(srcport) && error("$srcnode has no free ports")
    isnothing(dstport) && error("$dstnode has no free ports")
    isconnected(srcnode, srcport) && error("Port '$srcport' in $srcnode is already connected")
    isconnected(dstnode, dstport) && error("Port '$dstport' in $dstnode is already connected")

    return g.add_edge!(bg, (srcnode, srcport), (dstnode, dstport))
end

"""
    disconnect!(bg::BondGraph, node1, node2)

Remove the bond connecting `node1` and `node2`. The order of nodes does not matter.
"""
function disconnect!(bg::BondGraph, node1::AbstractNode, node2::AbstractNode)
    # rem_edge! removes the bond regardless of the direction of the bond
    return g.rem_edge!(bg, node1, node2)
end


"""
    swap!(bg::BondGraph, oldnode, newnode)

Remove `oldnode` from bond graph `bg` and replace it with `newnode`. The new node
will have the same connections (bonds) as the original model.

`newnode` must have a greater or equal number of ports as `oldnode`.
"""
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

"""
    insert_node!(bg::BondGraph, bond, newnode)
    insert_node!(bg::BondGraph, (node1, node2), newnode)

Inserts `newnode` between two existing connected nodes. The direction of the original bond
is preserved.

Supply either the two nodes as a tuple, or the bond that connects them in `bg`.
"""
function insert_node!(bg::BondGraph, bond::Bond, newnode::AbstractNode)
    src = srcnode(bond)
    dst = dstnode(bond)

    disconnect!(bg, src, dst)

    try
        add_node!(bg, newnode)
        connect!(bg, newnode, dst) # dst first for TF component
        connect!(bg, src, newnode)
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

"""
    merge_nodes!(bg::BondGraph, node1, node2; junction=EqualEffort())

Combine two copies of the same component in `bg` by adding a `junction` and connecting the
neighbours of `node1` and `node2` to the new junction.

Merging nodes this way means there is only one component representing a system element, and
all other nodes connect to the component via the new junction.
"""
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

"""
    simplify_junctions!(bg::BondGraph; remove_redundant=true, squash_identical=true)

Remove unnecessary or redundant Junctions from bond graph `bg`.

If `remove_redundant` is true, junctions that have zero or one neighbours are removed, and
junctions with two neighbours are squashed (connected components remain connected).

If `squash_identical` is true, connected junctions of the same type are squashed into a
single junction.
"""
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
