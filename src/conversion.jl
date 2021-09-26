# Convert a ReactionSystem type into a BondGraph type
function BondGraph(rs::ReactionSystem)
    bg = BondGraph(string(rs.name))

    all_reactions = reactions(rs)
    for (i, reaction) in enumerate(all_reactions)
        if i > 1 && is_reverse_off_previous(reaction, all_reactions[i-1])
            # If this reaction is an exact reverse of the previous reaction,
            # then together they form a bi-directional reaction pair.
            # Skip the second reaction
            continue
        end

        Re = Component(:Re, "R$i", numports=2)
        add_node!(bg, Re)
        half_equation!(bg, reaction.substrates, Re)
        half_equation!(bg, reaction.products, Re)
    end

    bg
end

is_reverse_off_previous(r1, r2) = Set(r1.substrates) == Set(r2.products) && Set(r2.substrates) == Set(r1.products)

function half_equation!(bg, species, Re)
    if length(species) > 1
        one_junction = Junction(:𝟏)
        add_node!(bg, one_junction)

        for spcs in species
            comp = Component(:Ce, stringify_species(spcs))
            add_node!(bg, comp)
            connect!(bg, comp, one_junction)
        end

        connect!(bg, one_junction, Re)
    else
        comp = Component(:Ce, stringify_species(species[1]))
        add_node!(bg, comp)
        connect!(bg, comp, Re)
    end
end

# removes "(t)" from the end of the species name
stringify_species(species) = string(species)[1:end-3]