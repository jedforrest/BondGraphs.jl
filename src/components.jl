@parameters t
const TIME = t
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
        e => "Flow"
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

# Chemical species (:Ce)
@parameters k R T
@variables q_1(t)
Ce_dict = Dict(
    :description => "Chemical species",
    :numports => 1,
    :parameters => OrderedDict(
        k => "Biochemical Constant; exp(mu_0/RT)/V",
        R => "Universal Gas Constant",
        T => "Temperature"
    ),
    :state_vars => OrderedDict(
        q => "Molar Quantity"
    ),
    :equations => [
        0 ~ R*T*log(k*q) - E[1],
        D(q) ~ F[1]
    ]
)

# Normalised chemical species (:ce)
@parameters k
@variables q(t)
ce_dict = Dict(
    :description => "Chemical species (normalised)",
    :numports => 1,
    :parameters => OrderedDict(
        k => "Biochemical Constant; exp(mu_0/RT)/V"
    ),
    :state_vars => OrderedDict(
        q => "Molar Quantity"
    ),
    :equations => [
        0 ~ log(k*q) - E[1],
        D(q) ~ F[1]
    ]
)

# Chemical reaction (:Re)
@parameters r R T
Re_dict = Dict(
    :description => "Biochemical reaction",
    :numports => 2,
    :parameters => OrderedDict(
        r => "Reaction Rate",
        R => "Universal Gas Constant",
        T => "Temperature"
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - r*(exp(E[1]/R/T) - exp(E[2]/R/T))
    ]
)

# Normalised chemical reaction (:re)
@parameters r R T
re_dict = Dict(
    :description => "Biochemical reaction (normalised)",
    :numports => 2,
    :parameters => OrderedDict(
        r => "Reaction Rate"
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - r*(exp(E[1]) - exp(E[2]))
    ]
)

standard_library = Dict(
    :R => R_dict,
    :C => C_dict,
    :I => I_dict,
    :Se => Se_dict,
    :Sf => Sf_dict,
    :TF => TF_dict,
    :Ce => Ce_dict,
    :ce => ce_dict,
    :Re => Re_dict,
    :re => re_dict
)