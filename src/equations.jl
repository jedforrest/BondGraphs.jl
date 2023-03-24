@variables t

# New simplification rules
#TODO: enable threaded simplifier option
#TODO: create tests for rule rewriting
exponent_rules = [
    @rule(exp(log(~x)) => ~x),
    @rule(log(exp(~x)) => ~x),
    @acrule(exp(~x + ~y) => exp(~x) * exp(~y)),
    @acrule(exp(~x * ~y) => exp(~y)^~x),
    @acrule(log(~x) + log(~y) => log(~x * ~y)),
    @acrule(log((~x)^(~a)) => ~a * log(~x)),
    @acrule(~a * exp(~b * log(~x)) => (~a) * (~x)^(~b)),
]
rw_exp = RestartedChain(exponent_rules)
rw_chain = RestartedChain([SymbolicUtils.default_simplifier(), rw_exp])
rewriter = Postwalk(rw_chain)

"""
    constitutive_relations(n::AbstractNode)

Return the constitutive relations (equations) for node `n`.

If `n` is a Junction, the flow and effort constraints are generated from its connections.
"""
constitutive_relations(n::AbstractNode) = equations(n)
function constitutive_relations(n::EqualEffort)
    if all(==(0), ports(n)) # all ports are empty
        return Equation[]
    end

    N = numports(n)
    @variables E(t)[1:N] F(t)[1:N]

    flow_constraint = [0 ~ sum(collect(F))]
    effort_constraints = [0 ~ E[1] - e for e in collect(E[2:end])]
    return vcat(flow_constraint, effort_constraints)
end
function constitutive_relations(n::EqualFlow)
    if all(==(0), ports(n)) # all ports are empty
        return Equation[]
    end

    N = numports(n)
    @variables E(t)[1:N] F(t)[1:N]

    W = weights(n)
    weighted_e = W .* collect(E)
    weighted_f = W .* collect(F)

    effort_constraint = [0 ~ sum(weighted_e)]
    flow_constraints = [0 ~ weighted_f[1] - f for f in collect(weighted_f[2:end])]
    return vcat(effort_constraint, flow_constraints)
end

"""
    constitutive_relations(bg::BondGraph; sub_defaults=false)

Generate the constitutive relations (equations) for bond graph `bg`.

The equations are symbolically derived from the equations of all the nodes and bonds in the
bond graph. If `sub_defaults` is true, the default parameter values for each component are
subbed into the equations.

NOTE: This creates a ModelingToolkit.ODESystem of the bond graph and returns only the
equations. If you want to numerically solve the bond graph equations, either use
[`simulate`](@ref) or create an ODESystem directly.
"""
function constitutive_relations(bg::BondGraph; sub_defaults=false)
    cr = simplify.(full_equations(ODESystem(bg)); expand=true, rewriter)
    # TODO: for some reason, must be simplified twice to work
    cr = simplify.(cr; rewriter=rewriter)
    if sub_defaults
        return _sub_defaults(cr, all_variables(bg))
    else
        return cr
    end
end
function constitutive_relations(bgn::BondGraphNode)
    return constitutive_relations(bgn.bondgraph)
end

@connector function MTKPort(; name)
    vars = @variables E(t) F(t) [connect = Flow]
    ODESystem(Equation[], t, vars, []; name=name)
end

function _sub_defaults(eqs, defaults)
    for (comp, var_dict) in defaults
        cname = name(comp)
        # This is sort of a hack to match BG vars to MTK vars
        # Does not include time-dependent variables
        function mtkvar(cname, var)
            name = Symbol("$(cname)₊$var")
            v, = @variables $name
            return v
        end
        namespaced_var_dict = Dict(mtkvar(cname, k) => v for (k, v) in var_dict)
        eqs = [substitute(eq, namespaced_var_dict) for eq in eqs]
    end
    simplify.(eqs; rewriter)
end

function get_connection_eq(b::Bond, subsystems)
    src, dst = b # (node, label)
    src_port = getproperty(subsystems[src[1]], Symbol("p$(src[2])"))
    dst_port = getproperty(subsystems[dst[1]], Symbol("p$(dst[2])"))
    connect(src_port, dst_port)
end

# AbstractNode
function ModelingToolkit.ODESystem(n::AbstractNode; name=name(n))
    N = numports(n)
    ps = [MTKPort(name=Symbol("p$i")) for i in 1:N]

    @variables E(t)[1:N] F(t)[1:N]
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

# BondGraph
function ModelingToolkit.ODESystem(m::BondGraph; simplify_eqs=true)
    # TODO check for disconnected ports/nodes
    (subsystems, connections) = get_subsys_and_connections(m)
    sys = compose_bg_model(subsystems, connections, m.name, simplify_eqs)
    structural_simplify(sys)
end

# BondGraphNode
function ModelingToolkit.ODESystem(bgn::BondGraphNode; name=name(bgn), simplify_eqs=false)
    N = numports(bgn)
    ps = [MTKPort(name=Symbol("p$i")) for i in 1:N]

    (subsystems, connections) = get_subsys_and_connections(bgn.bondgraph)
    es = [subsystems[comp].p1.E for comp in exposed(bgn)]
    fs = [subsystems[comp].p1.F for comp in exposed(bgn)]
    port_eqs = [
        [0 ~ p.E - E for (E, p) in zip(es, ps)]
        [0 ~ p.F + F for (F, p) in zip(fs, ps)]
    ]
    eqs = [connections; port_eqs]
    sys = compose_bg_model(subsystems, eqs, name, simplify_eqs)
    compose(sys, ps...)
end

function get_subsys_and_connections(bg::BondGraph)
    # Collect constitutive relations from components
    subsystems = OrderedDict(n => ODESystem(n) for n in nodes(bg))
    # Add constraints from bonds
    connections = [get_connection_eq(b, subsystems) for b in bonds(bg)]
    return subsystems, connections
end

function compose_bg_model(subsystems, eqs, name, simplify_eqs)
    @named _model = ODESystem(eqs, t; name=Symbol(name))
    model = compose(_model, collect(values(subsystems)))

    if simplify_eqs
        # NOTE: this breaks for DAEs
        # Something to do with missing variables from the variable map
        model = structural_simplify(model)
        neweqs = full_equations(model; simplify=true)
        neweqs = simplify.(neweqs; expand=true, rewriter=rewriter)
        @set! model.eqs = neweqs
        @set! model.substitutions = nothing
        return model
    else
        return model
    end
end


"""
    simulate(bg::BondGraph, tspan; u0=[], pmap=[], solver=Tsit5(), flag_ODE=true, kwargs...)

Simulate the bond graph model.

The keyword arguments are the same as for `ODEProblem` and `solve` in DifferentialEquations.
"""
function simulate(bg::BondGraph, tspan; u0=[], pmap=[], solver=Tsit5(), flag_ODE=true, kwargs...)
    # DAEs break custom model simplification, so skip this step in ODESystem
    sys = ODESystem(bg; simplify_eqs=flag_ODE)

    # If bg has control variables, need to allow union type for parameters
    use_union = has_controls(bg)

    # Check if problem is an ODE or DAE
    ODEProblemType = flag_ODE ? ODEProblem : ODAEProblem

    prob = ODEProblemType(sys, u0, tspan, pmap; use_union, kwargs...)
    return solve(prob, solver; kwargshandle=KeywordArgSilent)
end

# Custom post-processing of latex display for equations
# TODO: move inside cr functions
function Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation})
    ltx = latexify(x)

    # _+ becomes ₊
    ltx = replace(ltx, r"{\\_\+}" => s"_+")
    # dX₊q(t) becomes dqₓ(t)
    # ltx = replace(ltx, r"([^d\W]+){\\_\+}(\w+)" => s"\2_{\1}")
    # \mathrm{...} is removed
    # ltx = replace(ltx, r"\\mathrm{(.+?)}" => s"\1")

    print(io, "\$\$ " * ltx * " \$\$")
end
