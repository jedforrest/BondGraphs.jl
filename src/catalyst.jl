#TODO: include reaction rates in forming the bondgraph

# Convert a ReactionSystem type into a BondGraph type
function BondGraph(rs::ReactionSystem; chemostats=[])
    bg = BondGraph(rs.name)

    # Create disjoint reaction bondgraphs for each reaction in network
    all_reactions = reactions(rs)
    reaction_num = 1
    for (i, reaction) in enumerate(all_reactions)
        if i > 1 && _is_reverse_off_previous(reaction, all_reactions[i-1])
            # Skip the second reaction
            continue
        end

        Re = Component(:Re, Symbol("R$reaction_num"))
        add_node!(bg, Re)
        _half_equation!(bg, reaction.substrates, reaction.substoich, Re, chemostats)
        _half_equation!(bg, reaction.products, reaction.prodstoich, Re, chemostats)

        reaction_num += 1
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
_is_reverse_off_previous(r1, r2) = Set(r1.substrates) == Set(r2.products) && Set(r2.substrates) == Set(r1.products)

function _half_equation!(bg, species, stoich, Re, chemostats)
    species_names = _stringify_species.(species)

    if length(species) > 1
        one_junction = EqualFlow()
        add_node!(bg, one_junction)

        for (i, spcs) in enumerate(species_names)
            comp = spcs in chemostats ? Component(:Se, Symbol(spcs)) : Component(:Ce, Symbol(spcs))
            add_node!(bg, comp)
            connect!(bg, comp, one_junction)

            n = stoich[i]
            if n != 1
                tf = Component(:TF, "tf$(comp.name)"; n)
                insert_node!(bg, (comp, one_junction), tf)
            end
        end

        connect!(bg, one_junction, Re)
    else
        spcs = species_names[1]
        comp = spcs in chemostats ? Component(:Se, Symbol(spcs)) : Component(:Ce, Symbol(spcs))
        add_node!(bg, comp)
        connect!(bg, comp, Re)

        n = stoich[1]
        if n != 1
            tf = Component(:TF, "tf$(comp.name)"; n)
            insert_node!(bg, (comp, Re), tf)
        end
    end
end

# removes "(t)" from the end of the species name
_stringify_species(species) = string(species)[1:end-3]
