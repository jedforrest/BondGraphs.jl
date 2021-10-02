# New simplification rules
exponent_rules = [
    @rule(exp(log(~x)) => ~x),
    @acrule(exp(~x + ~y) => exp(~x) * exp(~y)),
    @acrule(exp(~x * ~y) => exp(~y)^~x),
]
rw_exp = RestartedChain(exponent_rules)
rw_chain = RestartedChain([SymbolicUtils.default_simplifier(),rw_exp])
rewriter = Postwalk(rw_chain)

# Equations
equations(n::Component) = n.equations
function equations(m::EqualEffort)
    if isempty(freeports(m))
        return Vector{Equation}([])
    end

    N = numports(m)
    @variables E[1:N](t) F[1:N](t)

    flow_constraint = [0 ~ sum(collect(F))]
    effort_constraints = [0 ~ E[1] - e for e in E[2:end]]
    return vcat(flow_constraint,effort_constraints)
end
function equations(m::EqualFlow)
    if isempty(freeports(m))
        return Vector{Equation}([])
    end
    N = numports(m)
    @variables E[1:N](t) F[1:N](t)

    weighted_e = m.weights.*collect(E)
    weighted_f = m.weights.*collect(F)

    effort_constraint = [0 ~ sum(weighted_e)]
    flow_constraints = [0 ~ weighted_f[1] - f for f in weighted_f[2:end]]
    return vcat(effort_constraint,flow_constraints)
end
function equations(m::BondGraph)
    return simplify.(ModelingToolkit.equations(de_system(m)),rewriter=rewriter)
end

function de_system(m::BondGraph)
    # Import equations from components, then substitute
    inv_x, inv_bs, inv_u = inverse_mappings(m)
    eqs = Vector{Equation}([])
    for c in components(m)
        new_eqs = equations(c)

        x_subs = Dict(x => inv_x[(c,x)] for x in state_vars(c))
        u_subs = Dict(p => inv_u[(c,p)] for p in params(c))

        N = numports(c)
        @variables E[1:N](t) F[1:N](t)
        e_subs = Dict(E[p.index] => e for (p,(e,_)) in inv_bs if p.node == c)
        f_subs = Dict(F[p.index] => f for (p,(_,f)) in inv_bs if p.node == c)
        
        sub_rules = merge(x_subs,u_subs,e_subs,f_subs)
        sub_eqs = [substitute(eq,sub_rules) for eq in new_eqs]
        append!(eqs,sub_eqs)
    end

    # Add constraints from bonds
    for (t,h) in bonds(m)
        (e_t,f_t) = inv_bs[t]
        (e_h,f_h) = inv_bs[h]
        push!(eqs, 0 ~ e_t - e_h)
        push!(eqs, 0 ~ f_t + f_h)
    end
    sys = ODESystem(eqs, t)
    return structural_simplify(sys)
end