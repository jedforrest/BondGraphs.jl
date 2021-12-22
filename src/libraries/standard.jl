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
    :equations => [0 ~ E[1]/R - F[1]]
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
    :state_vars => OrderedDict(
        q => "Generalised Position"
    ),
    :equations => [
        0 ~ q/C - E[1],
        D(q) ~ F[1]
    ]
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
    :state_vars => OrderedDict(
        p => "Generalised Momentum"
    ),
    :equations => [
        0 ~ p/L - F[1],
        D(p) ~ E[1]
    ]
)

# Source of effort (:Se)
@parameters e
Se_dict = Dict(
    :description => "Effort source",
    :numports => 1,
    :parameters => OrderedDict(
        e => "Effort"
    ),
    :equations => [0 ~ e - E[1]]
)

# Source of flow (:Sf)
@parameters f
Sf_dict = Dict(
    :description => "Flow source",
    :numports => 1,
    :parameters => OrderedDict(
        f => "Flow"
    ),
    :equations => [0 ~ f - F[1]]
)

# Transformer (:TF)
@parameters r
TF_dict = Dict(
    :description => "Linear Transformer",
    :numports => 2,
    :parameters => OrderedDict(
        r => "Winding ratio"
    ),
    :equations => [
        0 ~ E[2] - r*E[1],
        0 ~ F[1] - r*F[2]
    ]
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