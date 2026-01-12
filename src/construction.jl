"""
    add_comp!(bg::BondGraph, comps)

Add a comp to a bond graph `bg`.
"""
function add_comp!(bg::BondGraph, comp::AbstractElement)
    compname = name(comp)
    add_vertex!(bg.graph, compname, comp) || @warn "Component '$compname' already in model"
end

"""
    remove_comp!(bg::BondGraph, comps)

Remove a comp to a bond graph `bg`.
"""
function remove_comp!(bg::BondGraph, comp::AbstractElement)
    compname = name(comp)
    if !haskey(bg.graph, compname)
        @warn "Component '$compname' not in model"
        return false
    end
    # disconnect edge/bond before removing vertex
    for bond in filter(bond -> comp in bond, bonds(bg))
        disconnect!(bond)
    end
    vertex = code_for(bg.graph, compname)
    rem_vertex!(bg.graph, vertex)
end

"""
    connect!(bg::BondGraph, source_comp, destination_comp)
    connect!(bg::BondGraph, (source_comp, port_label), (destination_comp, port_label))

Connect two components together in the same bond graph. The bond direction is always from
`source_comp` to `destination_comp`. The port index of `source_comp` and `destination_comp`
can be optionally set.
"""
function connect!(bg::BondGraph, src::Union{Port,AbstractElement}, dst::Union{Port,AbstractElement})
    srcname = src isa Port ? parentname(src) : name(src)
    dstname = dst isa Port ? parentname(dst) : name(dst)
    haskey(bg.graph, srcname, dstname) && return true
    srcname in labels(bg.graph) || error("Component '$srcname' not found in bond graph")
    dstname in labels(bg.graph) || error("Component '$dstname' not found in bond graph")
    bond = Bond(src, dst)
    add_edge!(bg.graph, srcname, dstname, bond)
end

"""
    disconnect!(bg::BondGraph, comp1, comp2)

Remove the bond connecting `comp1` and `comp2`
"""
function disconnect!(bg::BondGraph, src::AbstractElement, dst::AbstractElement)
    srcname = name(src)
    dstname = name(dst)
    bond = bg[srcname, dstname]
    srcvertex = code_for(bg.graph, srcname)
    dstvertex = code_for(bg.graph, dstname)
    removed = rem_edge!(bg.graph, srcvertex, dstvertex)
    if removed
        disconnect!(bond)
    end
    removed
end

"""
    swap!(bg::BondGraph, oldcomp, newcomp)

Remove `oldcomp` from bond graph `bg` and replace it with `newcomp`. The new comp
will have the same connections (bonds) as the original model.

`newcomp` must have a greater or equal number of ports as `oldcomp`.
"""
function swap!(bg::BondGraph, oldcomp::AbstractElement, newcomp::AbstractElement)

    if numports(newcomp) < numports(oldcomp)
        @warn("New comp must have a greater or equal number of ports than the old comp")
        return false
    end

    # add new component to bond graph if not present already
    if !(newcomp in components(bg))
        add_comp!(bg, newcomp)
    end

    in_nbrs = inneighbor_comps(bg, oldcomp)
    out_nbrs = outneighbor_comps(bg, oldcomp)
    remove_comp!(bg, oldcomp)

    for i in in_nbrs
        in_comp = bg[i]
        connect!(bg, in_comp, newcomp)
    end
    for j in out_nbrs
        out_comp = bg[j]
        connect!(bg, newcomp, out_comp)
    end
    true
end

# """
#     insert_comp!(bg::BondGraph, bond, newcomp)
#     insert_comp!(bg::BondGraph, (comp1, comp2), newcomp)

# Inserts `newcomp` between two existing connected comps. The direction of the original bond
# is preserved.

# Supply either the two comps as a tuple, or the bond that connects them in `bg`.
# """
# function insert_comp!(bg::BondGraph, bond::Bond, newcomp::AbstractElement)
#     src = srccomp(bond)
#     dst = dstcomp(bond)

#     disconnect!(bg, src, dst)

#     try
#         add_comp!(bg, newcomp)
#         connect!(bg, src, newcomp)
#         connect!(bg, newcomp, dst)
#     catch e
#         # if connection fails, reconnect original bond
#         disconnect!(bg, src, newcomp)
#         disconnect!(bg, newcomp, dst)
#         connect!(bg, src, dst)
#         error(e)
#     end
# end
# function insert_comp!(bg::BondGraph, tuple::Tuple, newcomp::AbstractElement)
#     bonds = getbonds(bg, tuple)
#     isempty(bonds) && error("$(tuple[1]) and $(tuple[2]) are not connected")
#     insert_comp!(bg, bonds[1], newcomp)
# end

# """
#     merge_comps!(bg::BondGraph, comp1, comp2; junction=EqualEffort())

# Combine two copies of the same component in `bg` by adding a `junction` and connecting the
# neighbours of `comp1` and `comp2` to the new junction.

# Merging comps this way means there is only one component representing a system compent, and
# all other comps connect to the component via the new junction.
# """
# function merge_comps!(
#         bg::BondGraph,
#         comp1::Abstractcomp,
#         comp2::Abstractcomp;
#         junction = EqualEffort()
# )
#     comp1.type == comp2.type ||
#         error("$(comp1.name) must be the same type as $(comp2.name)")

#     # comp1 taken as the comp to keep
#     for nb in all_neighbors(bg, comp1)
#         junc = deepcopy(junction)
#         bond = getbonds(bg, comp1, nb)[1]
#         insert_comp!(bg, bond, junc)
#         swap!(bg, comp2, junc)
#     end
# end
# function merge_comps!(bg::BondGraph, comp1::Junction, comp2::Junction)
#     # comp1 taken as the comp to keep
#     # remove conflicting connections between junctions if they exist
#     disconnect!(bg, comp1, comp2)
#     shared_neighbors = intersect(all_neighbors(bg, comp1), all_neighbors(bg, comp2))
#     for shared_neighbor in shared_neighbors
#         disconnect!(bg, comp2, shared_neighbor)
#     end
#     swap!(bg, comp2, comp1)
# end

# """
#     simplify_junctions!(bg::BondGraph; remove_redundant=true, squash_identical=true)

# Remove unnecessary or redundant Junctions from bond graph `bg`.

# If `remove_redundant` is true, junctions that have zero or one neighbours are removed, and
# junctions with two neighbours are squashed (connected components remain connected).

# If `squash_identical` is true, connected junctions of the same type are squashed into a
# single junction.
# """
# function simplify_junctions!(
#         bg::BondGraph;
#         remove_redundant = true,
#         squash_identical = true
# )
#     junctions = filter(n -> n isa Junction, bg.comps)

#     # Removes junctions with 2 or less connected ports
#     if remove_redundant
#         for j in junctions
#             n_nbrs = length(all_neighbors(bg, j))
#             if n_nbrs == 2
#                 #srccomp = inneighbors(bg, j)[1]
#                 #dstcomp = outneighbors(bg, j)[1]
#                 comp1, comp2 = all_neighbors(bg, j)
#                 remove_comp!(bg, j)
#                 # bond direction may not be preserved here
#                 connect!(bg, comp1, comp2)
#             elseif n_nbrs < 2
#                 remove_comp!(bg, j)
#             end
#         end
#     end

#     # Squashes identical copies of the same junction type into one junction
#     if squash_identical
#         for j in junctions, nbr in all_neighbors(bg, j)

#             has_vertex(bg, j) || continue # in case j was removed
#             if type(j) == type(nbr)
#                 merge_comps!(bg, j, nbr)
#             end
#         end
#     end
#     bg
# end
