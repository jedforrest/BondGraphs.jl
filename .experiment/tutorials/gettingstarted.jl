# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

@variables e(t) f(t) p(t) q(t)
@parameters R C
phiR(e, f, R) = R * f - e
phiC(e, q, C) = q - C * e

r = DissipatorElement(phiR, [R], 1)
c = StaticStorageElement(phiC, [C], 1)
r(e, f)
c(e, f)

@named cap = Component(c)
@named res = Component(r)
@named j0 = Component(EqualEffort())

b1 = Bond(res, j0)
b2 = Bond(cap, j0)

bg = BondGraph("RC Circuit", [res, cap, j0], [b1, b2])
# bg = BondGraph("RC Circuit", [b1, b2]) # alternative

# using Graphs
# g = graph(bg)
# incidence_matrix(g)
# graphplot(bg)

system(res)
system(cap)
sys0 = system(j0)
equations(expand_connections(sys0))

sys = system(bg; simplify=false)
hierarchy(sys)

equations(sys)
equations(expand_connections(sys))
sys2 = structural_simplify(sys)  # will become mtkcompile

unknowns(sys2)
equations(sys2)
observed(sys2)
full_equations(sys2)

#############
@unpack q, C = sys4.capA
@unpack R = sys4.resA
sys4.capA.q
prob = ODEProblem(sys4, [sys4.capA.q => 1], (0., 10.), [sys4.resA.R => 2, sys4.capA.C => 1])
prob.ps
sol = solve(prob, Tsit5())
plot(sol)

#############
# TODO CONTINUE FROM HERE
C.C = 1
R.R = 2
constitutive_relations(model; sub_defaults=true)

tspan = (0., 10.)
u0 = [1] # initial value for C.q(t)
sol = simulate(model, tspan; u0)
plot(sol)

Is = Component(:Sf, "Is")
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
