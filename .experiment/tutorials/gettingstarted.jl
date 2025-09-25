# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

@variables e(t) f(t) p(t) q(t)
@parameters R C
phiR(e, f, R) = R * f - e
phiC(e, q, C) = C * q - e

r = DissipatorElement(phiR, [R], 1)
c = StaticStorageElement(phiC, [C], 1)
r(e, f)
c(e, f)

cap = Component(c, "C1")
res = Component(r, "R1")
j0 = Component(EqualEffort(), "j0")

b1 = Bond(res, j0)
b2 = Bond(j0, cap)

bg = BondGraph("RC Circuit", [res, cap, j0], [b1, b2])
# bg = BondGraph("RC Circuit", [b1, b2]) # alternative

# using Graphs
# g = graph(bg)
# incidence_matrix(g)
# graphplot(bg)

constitutive_relations(res)
constitutive_relations(cap)
constitutive_relations(j0)

system(res)
system(cap)

sys = system(bg; simplify=false)

equations(sys)
equations(expand_connections(sys))

# TODO CONTINUE FROM HERE
simplified_sys = structural_simplify(sys)  # will become mtkcompile


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
