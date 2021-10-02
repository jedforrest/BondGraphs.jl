function basis_vectors(m::BondGraph)
    tangent_space = [(x,D(x)) for x in state_vars(m)]
    # TODO: Fill in the varaibles below
    port_space = Dict()
    control_space = []
    return tangent_space, port_space, control_space
end

function bond_space(m::BondGraph)
    d = OrderedDict{Tuple{Num,Num},Port}()
    N = 2*length(bonds(m))
    @variables e[1:N](t) f[1:N](t)
    i = 0
    for b in bonds(m)
        i += 1
        d[(e[i], f[i])] = b.src
        i += 1
        d[(e[i], f[i])] = b.dst
    end
    return d
end
bond_space(m::AbstractNode) = OrderedDict{Tuple{Num,Num},Port}()


function control_space(m::AbstractNode)
    @parameters u[1:length(params(m))]
    return OrderedDict(u[i] => (m,p) for (i,p) in enumerate(params(m)))
end
function control_space(m::BondGraph)
    i = 0
    vars = Vector{Tuple{Model,Num}}([])
    for c in components(m), p in params(c)
        push!(vars,(c,p))
    end
    @parameters u[1:length(vars)]
    return OrderedDict(u[i] => v for (i,v) in enumerate(vars))
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
