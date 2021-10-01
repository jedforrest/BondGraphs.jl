function basis_vectors(m::BondGraph)
    tangent_space = [(x,D(x)) for x in state_vars(m)]
    # TODO: Fill in the varaibles below
    port_space = Dict()
    control_space = []
    return tangent_space, port_space, control_space
end

function bond_space(m::BondGraph)
    d = OrderedDict{Tuple{Num,Num},Port}()
    i = 0
    for b in bonds(m)
        i += 1
        d[(internal_effort(i),internal_flow(i))] = b.src
        i += 1
        d[(internal_effort(i),internal_flow(i))] = b.dst
    end
    return d
end
bond_space(m::AbstractNode) = OrderedDict{Tuple{Num,Num},Port}()


function control_space(m::AbstractNode)
    return OrderedDict(control_variable(i) => (m,p) for (i,p) in enumerate(params(m)))
end
function control_space(m::BondGraph)
    i = 0
    dict_states = OrderedDict{Num,Tuple{Model,Num}}()
    for c in components(m)
        for p in params(c)
            i += 1
            dict_states[control_variable(i)] = (c,p)
        end
    end
    return dict_states
end

function invert(d::OrderedDict)
    OrderedDict(v => k for (k,v) in d)
end

function inverse_mappings(m)
    return (
        invert(state_vars(m)),
        invert(bond_space(m)),
        invert(control_space(m))
    )
end
