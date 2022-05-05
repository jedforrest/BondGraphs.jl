@variables t

# New simplification rules
exponent_rules = [
    @rule(exp(log(~x)) => ~x),
    @rule(log(exp(~x)) => ~x),
    @acrule(exp(~x + ~y) => exp(~x) * exp(~y)),
    @acrule(exp(~x * ~y) => exp(~y)^~x),
    @acrule(log(~x) + log(~y) => log(~x * ~y)),
    # @rule(log((~x)^(~a)) => ~a * log(~x)),
    @rule(exp(~a * log(~x)) => (~x)^(~a))
]
rw_exp = RestartedChain(exponent_rules)
rw_chain = RestartedChain([SymbolicUtils.default_simplifier(), rw_exp])
rewriter = Postwalk(rw_chain)

# Constitutive relations
constitutive_relations(n::AbstractNode) = equations(n)
function constitutive_relations(n::EqualEffort)
    if all(freeports(n)) # all ports are empty
        return Equation[]
    end

    N = numports(n)
    @variables E[1:N](t) F[1:N](t)

    flow_constraint = [0 ~ sum(collect(F))]
    effort_constraints = [0 ~ E[1] - e for e in E[2:end]]
    return vcat(flow_constraint, effort_constraints)
end
function constitutive_relations(n::EqualFlow)
    if all(freeports(n)) # all ports are empty
        return Equation[]
    end

    N = numports(n)
    @variables E[1:N](t) F[1:N](t)

    W = weights(n)
    weighted_e = W .* collect(E)
    weighted_f = W .* collect(F)

    effort_constraint = [0 ~ sum(weighted_e)]
    flow_constraints = [0 ~ weighted_f[1] - f for f in weighted_f[2:end]]
    return vcat(effort_constraint, flow_constraints)
end
function constitutive_relations(m::BondGraph)
    return simplify.(full_equations(ODESystem(m)); expand = true, rewriter = rewriter)
end
function constitutive_relations(bgn::BondGraphNode)
    return constitutive_relations(bgn.bondgraph)
end

@connector function MTKPort(; name)
    vars = @variables E(t) F(t) [connect = Flow]
    ODESystem(Equation[], t, vars, []; name = name)
end

function get_connection_eq(b::Bond,subsystems)
    (s,d) = b
    src_port = getproperty(subsystems[s.node], Symbol("p$(s.index)"))
    dst_port = getproperty(subsystems[d.node], Symbol("p$(d.index)"))
    connect(src_port, dst_port)
end

function ModelingToolkit.ODESystem(n::AbstractNode; name=name(n))
    N = numports(n)
    ps = [MTKPort(name=Symbol("p$i")) for i in 1:N]

    @variables E[1:N](t) F[1:N](t)
    e_sub_rules = Dict(E[i] => ps[i].E for i in 1:N)
    f_sub_rules = Dict(F[i] => ps[i].F for i in 1:N)
    u_sub_rules = Dict(u => controls(n)[u](t) for u in keys(controls(n)))

    sub_rules = merge(e_sub_rules, f_sub_rules, u_sub_rules)
    eqs = Equation[substitute(eq, sub_rules) for eq in constitutive_relations(n)]

    # mtk vars are stored as keys
    _params = collect(keys(parameters(n)))
    _ctrls = collect(keys(controls(n)))
    _globals = collect(keys(globals(n)))
    _states = collect(keys(states(n)))

    sys = ODESystem(eqs, t, _states, [_params; _ctrls; _globals];
        name=Symbol(name), defaults=all_variables(n), controls=_ctrls)
    return compose(sys, ps...)
end
function ModelingToolkit.ODESystem(m::BondGraph; simplify_eqs = true)
    (subsystems, connections) = get_subsys_and_connections(m)
    compose_bg_model(subsystems, connections, m.name, simplify_eqs)
end
function ModelingToolkit.ODESystem(bgn::BondGraphNode; name=name(bgn), simplify_eqs = false)
    N = numports(bgn)
    ps = [MTKPort(name = Symbol("p$i")) for i in 1:N]

    (subsystems, connections) = get_subsys_and_connections(bgn.bondgraph)
    es = [subsystems[c].p1.E for c in exposed(bgn)]
    fs = [subsystems[c].p1.F for c in exposed(bgn)]
    port_eqs = [
        [0 ~ p.E - E for (E,p) in zip(es,ps)];
        [0 ~ p.F + F for (F,p) in zip(fs,ps)]
    ]
    eqs = [connections; port_eqs]
    sys = compose_bg_model(subsystems, eqs, name, simplify_eqs)
    compose(sys, ps...)
end

function get_subsys_and_connections(bg::BondGraph)
    uniquenames = create_unique_names(bg.nodes)
    # Collect constitutive relations from components
    subsystems = OrderedDict(n => ODESystem(n; name=uniquenames[n]) for n in nodes(bg))
    # Add constraints from bonds
    connections = [get_connection_eq(b, subsystems) for b in bonds(bg)]
    return subsystems, connections
end

function compose_bg_model(subsystems, eqs, name, simplify_eqs)
    @named _model = ODESystem(eqs, t; name = Symbol(name))
    model = compose(_model, collect(values(subsystems)))

    if simplify_eqs
        return structural_simplify(model)
    else
        return model
    end
end

function create_unique_names(nodes)
    nodenames = name.(nodes)
    newnames = Dict()
    for n in nodes
        _name = name(n)
        common_named_nodes = nodes[findall(==(_name), nodenames)]
        if length(common_named_nodes) > 1
            for (i, _n) in enumerate(common_named_nodes)
                newnames[_n] = _n.name * "$i"
            end
        else
            newnames[n] = n.name
        end
    end
    newnames
end

function simulate(m::BondGraph, tspan; u0 = [], pmap = [], probtype::Symbol = :default, kwargs...)
    sys = ODESystem(m)
    flag_ODE = !any([isequal(eq.lhs, 0) for eq in ModelingToolkit.equations(sys)])
    if probtype == :ODE || flag_ODE
        prob = ODEProblem(sys, u0, tspan, pmap)
    else
        prob = ODAEProblem(sys, u0, tspan, pmap)
    end
    return solve(prob; kwargs...)
end

# Custom post-processing of latex display for equations
# TODO: may remove
function Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation})
    ltx = latexify(x)

    # _+ becomes ₊
    ltx = replace(ltx, r"{\\_\+}" => s"_+")
    # dX₊q(t) becomes dqₓ(t)
    # ltx = replace(ltx, r"([^d\W]+){\\_\+}(\w+)" => s"\2_{\1}")
    # \mathrm{...} is removed
    ltx = replace(ltx, r"\\mathrm{(.+?)}" => s"\1")

    print(io, ltx)
end
