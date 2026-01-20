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
    !haskey(bg.graph, srcname, dstname) && return false
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

"""
    insert_comp!(bg::BondGraph, bond, newcomp)
    insert_comp!(bg::BondGraph, (comp1, comp2), newcomp)

Inserts `newcomp` between two existing connected comps. The direction of the original bond
is preserved.

Supply either the two comps as a tuple, or the bond that connects them in `bg`.
"""
function insert_comp!(bg::BondGraph, bond::Bond, newcomp::AbstractElement)
    srcname, dstname = componentnames(bond)
    src = bg[srcname]
    dst = bg[dstname]

    disconnect!(bg, src, dst)

    try
        add_comp!(bg, newcomp)
        connect!(bg, src, newcomp)
        connect!(bg, newcomp, dst)
    catch e
        # if connection fails, reconnect original bond
        disconnect!(bg, src, newcomp)
        disconnect!(bg, newcomp, dst)
        connect!(bg, src, dst)
        error(e)
    end
end
function insert_comp!(bg::BondGraph, tuple::Tuple{AbstractElement, AbstractElement}, newcomp::AbstractElement)
    srcname, dstname = name.(tuple)
    if !haskey(bg.graph, srcname, dstname)
        error("'$(srcname)' and '$(dstname)' are not connected")
    end
    insert_comp!(bg, bg[srcname, dstname], newcomp)
end

# TODO extend for ParametricJunction
function merge_junctions!(bg::BondGraph, comp1::T, comp2::T) where {T <: NonParametricJunction}
    # comp1 taken as the comp to keep
    # remove conflicting connections between junctions if they exist
    disconnect!(bg, comp1, comp2)
    all_nbrs1 = all_neighbor_labels(bg.graph, name(comp1))
    all_nbrs2 = all_neighbor_labels(bg.graph, name(comp2))
    shared_neighbors = intersect(all_nbrs1, all_nbrs2)
    for shared_neighbor in shared_neighbors
        shared_nbr_comp = bg[shared_neighbor]
        disconnect!(bg, comp2, shared_nbr_comp)
    end
    swap!(bg, comp2, comp1)
end

"""
    simplify_junctions!(bg::BondGraph; remove_redundant=true, squash_identical=true)

Remove unnecessary or redundant Junctions from bond graph `bg`.

If `remove_redundant` is true, junctions that have zero or one neighbours are removed, and
junctions with two neighbours are squashed (connected components remain connected).

If `squash_identical` is true, connected junctions of the same type are squashed into a
single junction.
"""
function simplify_junctions!(
        bg::BondGraph;
        remove_redundant = true,
        squash_identical = true
)
    junctions = filter(x -> x isa EqualEffort || x isa EqualFlow, components(bg))

    # Removes junctions with 2 or less connected ports
    if remove_redundant
        for j in junctions
            all_nbrs = all_neighbor_labels(bg.graph, name(j))
            if length(all_nbrs) == 2
                comp1, comp2 = all_nbrs
                remove_comp!(bg, j)
                # bond direction may not be preserved here
                connect!(bg, bg[comp1], bg[comp2])
            elseif length(all_nbrs) < 2
                remove_comp!(bg, j)
            end
        end
    end

    # Squashes identical copies of the same junction type into one junction
    # FIXME
    if squash_identical
        for j in junctions
            all_nbrs = collect(all_neighbor_labels(bg.graph, name(j)))
            for nbr in all_nbrs
                g = bg.graph
                has_vertex(g, code_for(g, name(j))) || continue # in case j was removed
                nbrcomp = bg[nbr]
                if typeof(j) == typeof(nbrcomp)
                    merge_junctions!(bg, j, nbrcomp)
                end
            end
        end
    end

    # TODO merge TF/GY components

    bg
end
