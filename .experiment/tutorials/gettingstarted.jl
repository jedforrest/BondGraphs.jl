# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

model = BondGraph(name="RC Circuit")

C = StorageElement(:C)

R = Component(:R)
kvl = EqualEffort()

add_node!(model, [C, R, kvl])
connect!(model, R, kvl)
connect!(model, C, kvl)
model

using Graphs
incidence_matrix(model)

constitutive_relations(model)

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
