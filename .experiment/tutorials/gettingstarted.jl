# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

@variables e f p q R C
@variables ee[1:3] ff[1:3]
phiR(e, f, R) = R * f - e
phiC(e, q, C) = C * q - e

r = DissipatorElement(phiR, [R], 1)
c = StaticStorageElement(phiC, [C], 1)
z = EqualEffort()
r(e, f)
c(e, f)
z(collect(ee), collect(ff))

cap = Component(c, "C1")
res = Component(r, "R1")
j0 = Component(z, "j0")

subtype(cap)
parameters(cap)

b1 = Bond(res, j0)
b2 = Bond(j0, cap)

bg = BondGraph("RC Circuit", [res, cap, j0], [b1, b2])
# bg = BondGraph("RC Circuit", [b1, b2]) # alternative

using Graphs
g = graph(bg)
incidence_matrix(g)

constitutive_relations(res)
constitutive_relations(cap)
constitutive_relations(j0)

equations(bg)



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
