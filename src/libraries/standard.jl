module StandardLibrary

using ModelingToolkit
using OrderedCollections

export standard_library

@parameters t
D = Differential(t)

@variables E[1:2](t) F[1:2](t)

# Linear resistance (:R)
@parameters R
R_dict = Dict(
    :description => "Generalised Linear Resistor",
    :numports => 1,
    :parameters => OrderedDict(
        R => "Resistance"
    ),
    :equations => [0 ~ E[1] - R*F[1]],
    :defaults => Dict(R => 1.)
)

# Linear capacitor (:C)
@parameters C
@variables q(t)
C_dict = Dict(
    :description => "Generalised Linear Resistor",
    :numports => 1,
    :parameters => OrderedDict(
        C => "Capacitance"
    ),
    :states => OrderedDict(
        q => "Generalised Position"
    ),
    :equations => [
        0 ~ q/C - E[1],
        D(q) ~ F[1]
    ],
    :defaults => Dict(C => 1., q => 0.)
)

# Linear inductance (:I)
@parameters L
@variables p(t)
I_dict = Dict(
    :description => "Generalised Linear Inductor",
    :numports => 1,
    :parameters => OrderedDict(
        L => "Inductance"
    ),
    :states => OrderedDict(
        p => "Generalised Momentum"
    ),
    :equations => [
        0 ~ p/L - F[1],
        D(p) ~ E[1]
    ],
    :defaults => Dict(L => 1., p => 0.)
)

# Source of effort (:Se)
@parameters e
Se_dict = Dict(
    :description => "Effort source",
    :numports => 1,
    :parameters => OrderedDict(
        e => "Effort"
    ),
    :equations => [0 ~ e - E[1]],
    :defaults => Dict(e => 1.)
)

# Source of flow (:Sf)
@parameters f
Sf_dict = Dict(
    :description => "Flow source",
    :numports => 1,
    :parameters => OrderedDict(
        f => "Flow"
    ),
    :equations => [0 ~ f - F[1]],
    :defaults => Dict(f => 1.)
)

# Transformer (:TF)
@parameters n
TF_dict = Dict(
    :description => "Linear Transformer",
    :numports => 2,
    :parameters => OrderedDict(
        n => "Ratio"
    ),
    :equations => [
        0 ~ E[2] - n * E[1],
        0 ~ F[1] - n * F[2]
    ],
    :defaults => Dict(n => 1.)
)

const standard_library = Dict(
    :R => R_dict,
    :C => C_dict,
    :I => I_dict,
    :Se => Se_dict,
    :Sf => Sf_dict,
    :TF => TF_dict
)

end