# Examples

Interactive Jupyter Notebook versions of these tutorials can be found on [GitHub](https://github.com/jedforrest/BondGraphsTutorials).

## Simple Electric Circuit

This is a reduced copy-and-paste version of the electric circuit tutorial in the [Getting Started](@ref) section.

![](assets/simple_electric_circuit_current.png)


```@example
using BondGraphs
#using ModelingToolkit: @register_symbolic
using Plots

model = BondGraph("RC Circuit")
C = Component(:C)
R = Component(:R)
Is = Component(:Sf, "Is")
kvl = EqualEffort()

add_node!(model, [C, R, Is, kvl])
connect!(model, R, kvl)
connect!(model, C, kvl)
connect!(model, Is, kvl)

C.C = 1
R.R = 2
constitutive_relations(model; sub_defaults=true) |> display

u0 = [1]
p = plot()
for i in 1:4
    Is.fs = t -> cos(i * t)
    sol = simulate(model, (0., 5.); u0)
    plot!(p, sol, label = "f(t) = cos($(i)t)", lw=2)
end
plot(p)
```

## Biochemical Reaction Networks

## SERCA Pump

## Electrochemical system - Ion Transport

## Custom Components - Enzyme Catalysed Reactions

