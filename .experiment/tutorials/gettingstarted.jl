# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

include("../standardlibrary.jl")
using .Library

@named res = Library.Resistor()
@named cap = Library.Capacitor()

# R: e <=> f
r = DissipatorElement(res)
# C: e <=> q
c = StaticStorageElement(cap)

@named rcomp = Component(r)
@named ccomp = Component(c)
@named j0 = Component(EqualEffort())

b1 = Bond(rcomp, j0)
b2 = Bond(ccomp, j0)

bg = BondGraph("RC Circuit", [rcomp, ccomp, j0], [b1, b2])

system(rcomp)
system(ccomp)
sys0 = system(j0)
equations(expand_connections(sys0))

j0.element
j0.ports

sys0

basesys = j0.element.sys
portsys = system.(j0.ports, namespaced=false)
portsys[1]
sys = compose(basesys, portsys)

# TODO cleanup this code (its not very clear)
e, f = getefforts(sys), getflows(sys)
inner_connection_eqs = junc.element(e, f)
extend(System(inner_connection_eqs, t; name=sys.name), sys)

sys0 = system(j0)
equations(expand_connections(sys0))

b = b1
srcport_sys = system(b.src)
dstport_sys = system(b.dst)
ModelingToolkit.connect(srcport_sys, dstport_sys)

bgsys = system(bg; simplify=false)
hierarchy(bgsys)

###
comps = components(bg)
subsyss = system.(comps)

conn_eqns = connection_equation.(bg.bonds)
basesys = ODESystem(conn_eqns, t, name=bg.name)

for sub in subsyss
    eq = equations(expand_connections(sub))
    println.(eq)
    println()
end

nameof.(subsyss)
conn_eqns

equations(expand_connections(basesys))

sys = System(conn_eqns, t, name=bg.name; systems=subsyss)
equations(sys)
equations(expand_connections(sys))

sys2 = structural_simplify(sys)
###

equations(sys2)

equations(bgsys)
equations(expand_connections(bgsys))
sys2 = structural_simplify(bgsys)  # will become mtkcompile in later update

unknowns(sys2)
equations(sys2)
observed(sys2)
full_equations(sys2)

#############
prob = ODEProblem(sys2, [sys2.cap.q => 1], (0., 10.), [sys2.res.R => 2, sys2.cap.C => 1])
prob.ps
sol = solve(prob, Tsit5())
plot(sol)

#############
# using Graphs
# g = graph(bg)
# incidence_matrix(g)
# graphplot(bg)
#############

# TODO CONTINUE FROM HERE

Is = Component(FlowSource(), name=:Is)
add_node!(model, Is)
connect!(model, Is, kvl)
plot(model)

Is.fs = t -> sin(2t)
constitutive_relations(model; sub_defaults=true)

sol = simulate(model, tspan; u0)
plot(sol)

using ModelingToolkit
@register_symbolic f(t)
Is.fs = t -> f(t)

f(t) = t % 2 <= 1 ? 0 : 1 # repeating square wave
sol = simulate(model, tspan; u0)
plot(sol)
