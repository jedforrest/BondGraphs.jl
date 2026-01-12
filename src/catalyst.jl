#TODO: include reaction rates in forming the bondgraph

"""
    merge_comps!(bg::BondGraph, comp1, comp2; junction=EqualEffort())

Combine two copies of the same component in `bg` by adding a `junction` and connecting the
neighbours of `comp1` and `comp2` to the new junction.

Merging comps this way means there is only one component representing a system component, and
all other comps connect to the component via the new junction.
"""
# FIXME this is mainly useful for catalyst/chemical bond graphs
function merge_comps!(bg::BondGraph, comp1::T, comp2::T; junction = EqualEffort()) where {T <: AbstractElement}
    comp1.type == comp2.type || error("Components must have the same type")

    # comp1 taken as the comp to keep
    for nb in all_neighbors(bg, comp1)
        junc = deepcopy(junction)
        bond = getbonds(bg, comp1, nb)[1]
        insert_comp!(bg, bond, junc)
        swap!(bg, comp2, junc)
    end
end

"""
    BondGraph(rs::ReactionSystem; chemostats=[])

Convert a Catalyst.ReactionSystem into a BondGraph.

`chemostats` are chemical species with fixed concentrations. In bond graph terms, these are
"SCe" types (chemical energy sources) instead of "Ce" types (chemical energy store).
"""
function BondGraph(rs::ReactionSystem; chemostats = [])
    bg = BondGraph(nameof(rs))

    re_num = Ref(1)
    tf_num = Ref(1)

    # Create disjoint reaction bondgraphs for each reaction in network
    all_reactions = reactions(rs)
    for (i, reaction) in enumerate(all_reactions)
        if i > 1 && _is_reverse_off_previous(reaction, all_reactions[i - 1])
            # Skip the second reaction
            continue
        end

        Re = Component(:Re, Symbol("R$(re_num[])"))
        add_node!(bg, Re)
        _half_equation!(
            bg,
            reaction.substrates,
            reaction.substoich,
            Re,
            chemostats,
            tf_num,
            port = (Re, 1)
        )
        _half_equation!(
            bg,
            reaction.products,
            reaction.prodstoich,
            Re,
            chemostats,
            tf_num,
            port = (Re, 2)
        )

        re_num[] += 1
    end

    # Combine common species across reactions
    species_names = _stringify_species.(species(rs))
    for spcs_name in species_names
        spcs_nodes = getproperty(bg, Symbol(spcs_name))
        spcs_nodes isa Vector || continue
        for node in spcs_nodes[2:end]
            merge_nodes!(bg, spcs_nodes[1], node)
        end
    end

    simplify_junctions!(bg)
end

# If this reaction is an exact reverse of the previous reaction,
# then together they form a bi-directional reaction pair.
function _is_reverse_off_previous(r1, r2)
    Set(r1.substrates) == Set(r2.products) && Set(r2.substrates) == Set(r1.products)
end

function _half_equation!(bg, species, stoich, Re, chemostats, tf_num; port = Re)
    species_names = _stringify_species.(species)

    if length(species) > 1
        one_junction = EqualFlow()
        add_node!(bg, one_junction)

        for (i, spcs) in enumerate(species_names)
            comp = spcs in chemostats ? Component(:SCe, Symbol(spcs)) :
                   Component(:Ce, Symbol(spcs))
            add_node!(bg, comp)
            connect!(bg, comp, one_junction)

            n = stoich[i]
            n != 1 && _insert_tf!(bg, comp, one_junction, n, tf_num)
        end

        connect!(bg, one_junction, port)
    else
        spcs = species_names[1]
        comp = spcs in chemostats ? Component(:SCe, Symbol(spcs)) :
               Component(:Ce, Symbol(spcs))
        add_node!(bg, comp)
        connect!(bg, comp, port)

        n = stoich[1]
        n != 1 && _insert_tf!(bg, comp, Re, n, tf_num)
    end
end

# removes "(t)" from the end of the species name
_stringify_species(species) = string(species)[1:(end - 3)]

function _insert_tf!(bg, node1, node2, n, tf_num)
    tf = Component(:TF, "tf$(tf_num[])"; n)
    insert_node!(bg, (node1, node2), tf)
    tf_num[] += 1
end
