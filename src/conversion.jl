# Convert a ReactionSystem type into a BondGraph type
function BondGraph(rs::ReactionSystem; chemostats=[])
    bg = BondGraph(string(rs.name))

    # Create disjoint reaction bondgraphs for each reaction in network
    all_reactions = reactions(rs)
    reaction_num = 1
    for (i, reaction) in enumerate(all_reactions)
        if i > 1 && is_reverse_off_previous(reaction, all_reactions[i-1])
            # Skip the second reaction
            continue
        end

        Re = Component(:Re, "R$reaction_num", numports=2)
        add_node!(bg, Re)
        half_equation!(bg, reaction.substrates, reaction.substoich, Re, chemostats)
        half_equation!(bg, reaction.products, reaction.prodstoich, Re, chemostats)

        reaction_num += 1
    end

    # Combine common species across reactions
    species_names = stringify_species.(species(rs))
    for spcs_name in species_names
        spcs_nodes = getnodes(bg, spcs_name)
        for node in spcs_nodes[2:end]
            merge_nodes!(bg, spcs_nodes[1], node)
        end
    end

    bg
end

# If this reaction is an exact reverse of the previous reaction,
# then together they form a bi-directional reaction pair.
is_reverse_off_previous(r1, r2) = Set(r1.substrates) == Set(r2.products) && Set(r2.substrates) == Set(r1.products)

function half_equation!(bg, species, stoich, Re, chemostats)
    species_names = stringify_species.(species)

    if length(species) > 1
        one_junction = Junction(:𝟏)
        add_node!(bg, one_junction)

        for (i, spcs) in enumerate(species_names)
            comp = spcs in chemostats ? Component(:Se, spcs) : Component(:Ce, spcs)
            add_node!(bg, comp)
            connect!(bg, comp, one_junction)

            if stoich[i] != 1
                tf = Component(:TF, "$(stoich[i])", numports=2)
                insert_node!(bg, (comp, one_junction), tf)
            end
        end

        connect!(bg, one_junction, Re)
    else
        spcs = species_names[1]
        comp = spcs in chemostats ? Component(:Se, spcs) : Component(:Ce, spcs)
        add_node!(bg, comp)
        connect!(bg, comp, Re)

        if stoich[1] != 1
            tf = Component(:TF, "$(stoich[1])", numports=2)
            insert_node!(bg, (comp, Re), tf)
        end
    end
end

# removes "(t)" from the end of the species name
stringify_species(species) = string(species)[1:end-3]