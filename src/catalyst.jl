#TODO: include reaction rates in forming the bondgraph

"""
    merge_comps!(bg::BondGraph, comp1, comp2; junction=EqualEffort)

Combine two copies of the same component in `bg` by adding a `junction` and connecting the
neighbours of `comp1` and `comp2` to the new junction.

Merging comps this way means there is only one component representing a system component, and
all other comps connect to the component via the new junction.
"""
function merge_comps!(bg::BondGraph, comp1::T, comp2::T;
        junctiontype::Type{<:JunctionStructure} = EqualEffort) where {T <: BondElement}

    # comp1 taken as the comp to keep
    comp1_name = name(comp1)
    in_nbrs = inneighbor_comps(bg, comp1)
    out_nbrs = outneighbor_comps(bg, comp1)

    nbr_bonds = Bond[]
    for i in in_nbrs
        push!(nbr_bonds, bg[i, comp1_name])
    end
    for o in out_nbrs
        push!(nbr_bonds, bg[comp1_name, o])
    end

    juncname = "$(name(comp1))_$junctiontype"
    junc = junctiontype(; name=Symbol(juncname))

    for bond in nbr_bonds
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
function BondGraph(rn::ReactionSystem; chemostats = [], simplify=true)
    bg = BondGraph(; name=nameof(rn))

    # counter for added components
    comp_counter = counter([:re, :tf])

    # add species and zero junctions
    species_names = tosymbol.(species(rn), escape=false)
    for species_name in species_names
        species_comp = if species_name in chemostats
            chemostat(; name=species_name)
        else
            chemicalspecies(; name=species_name)
        end
        add_comp!(bg, species_comp)

        # zero junction for this species
        zero_junc_name = Symbol("0_$species_name")
        zero_junc = EqualEffort(; name=zero_junc_name)
        add_comp!(bg, zero_junc)

        connect!(bg, species_comp, zero_junc)
    end

    # Create disjoint reaction components for each reaction in network
    all_reactions = reactions(rn)
    for (i, rxn) in enumerate(all_reactions)
        if i > 1 && _is_reverse_off_previous(rxn, all_reactions[i - 1])
            continue  # Skip the second reaction
        end

        # Reaction component
        re_name = Symbol("Re$(comp_counter[:re])")
        Re = reaction(; name=re_name)
        add_comp!(bg, Re)
        inc!(comp_counter, :re)

        # One junctions for each reaction side
        one_junc_subt = EqualFlow(; name=Symbol("1_$(re_name)s"))
        one_junc_prod = EqualFlow(; name=Symbol("1_$(re_name)p"))
        add_comp!(bg, one_junc_subt)
        add_comp!(bg, one_junc_prod)

        # connect Reaction to One junctions
        connect!(bg, one_junc_subt, Re)
        connect!(bg, Re, one_junc_prod)

        # add species for each side of the reaction
        for j = 1:2
            if j == 1
                sps = rxn.substrates
                stoich = rxn.substoich
                one_junc = one_junc_subt
            else
                sps = rxn.products
                stoich = rxn.prodstoich
                one_junc = one_junc_prod
            end

            # connect to zero junctions already in the model
            species_names = tosymbol.(sps, escape=false)
            for (k, species_name) in enumerate(species_names)
                zero_junc = bg[Symbol("0_$species_name")]
                connect!(bg, zero_junc, one_junc)

                n = stoich[k]
                if n != 1
                    _insert_tf!(bg, zero_junc, one_junc, n, comp_counter[:tf])
                    inc!(comp_counter, :tf)
                end
            end
        end

    end

    # Combine common species across reactions
    species_names = tosymbol.(species(rn), escape=false)
    for spcs_name in species_names
        spcs_nodes = bg[spcs_name]
        if spcs_nodes isa Vector
            for node in spcs_nodes[2:end]
                merge_comps!(bg, spcs_nodes[1], node)
            end
        end
    end

    # Remove unnecessary junctions i.e. One junctions with two ports
    if simplify
        simplify_junctions!(bg, squash_identical=false)
    end

    bg
end

# If this reaction is an exact reverse of the previous reaction,
# then together they form a bi-directional reaction pair.
function _is_reverse_off_previous(r1, r2)
    Set(r1.substrates) == Set(r2.products) && Set(r2.substrates) == Set(r1.products)
end

function _insert_tf!(bg, node1, node2, n, tf_num)
    tfname = Symbol("TF$tf_num")
    tf = stoichiometry(n, name=tfname)
    insert_comp!(bg, (node1, node2), tf)
end
