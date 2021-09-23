@parameters t
D = Differential(t)

# Linear resistance (:R)
@parameters r
@variables E_1(t) F_1(t)
R_dict = Dict(
    :description => "Generalised Linear Resistor",
    :numports => 1,
    :params => Dict(
        r => "Resistance"
    ),
    :equations => [0 ~ E_1/r - F_1]
)

# Linear capacitor (:C)
@parameters C
@variables E_1(t) F_1(t) q_1(t)
C_dict = Dict(
    :description => "Generalised Linear Resistor",
    :numports => 1,
    :params => Dict(
        C => "Capacitance"
    ),
    :equations => [
        0 ~ q_1/C - E_1,
        D(q_1) ~ F_1
    ]
)

# Linear inductance (:I)
@parameters L
@variables E_1(t) F_1(t) p_1(t)
I_dict = Dict(
    :description => "Generalised Linear Inductor",
    :numports => 1,
    :params => Dict(
        L => "Inductance"
    ),
    :equations => [
        0 ~ p_1/L - F_1,
        D(p_1) ~ E_1
    ]
)

# Source of effort (:Se)
@parameters e
@variables E_1(t)
Se_dict = Dict(
    :description => "Effort source",
    :numports => 1,
    :params => Dict(
        e => "Effort"
    ),
    :equations => [0 ~ e - E_1]
)

# Source of flow (:Sf)
@parameters f
@variables F_1(t)
Sf_dict = Dict(
    :description => "Flow source",
    :numports => 1,
    :params => Dict(
        e => "Flow"
    ),
    :equations => [0 ~ f - F_1]
)

# Transformer (:TF)
@parameters r
@variables E_1(t) F_1(t) E_2(t) F_2(t)
TF_dict = Dict(
    :description => "Linear Transformer",
    :numports => 2,
    :params => Dict(
        r => "Winding ratio"
    ),
    :equations => [
        0 ~ E_2 - r*E_1,
        0 ~ F_1 - r*F_2
    ]
)

# Chemical species (:Ce)
@parameters k R T
@variables E_1(t) F_1(t) q_1(t)
Ce_dict = Dict(
    :description => "Chemical species",
    :numports => 1,
    :params => Dict(
        k => "Biochemical Constant; exp(mu_0/RT)/V",
        R => "Universal Gas Constant",
        T => "Temperature"
    ),
    :equations => [
        0 ~ R*T*log(k*q_1) - E_1,
        D(q_1) ~ F_1
    ]
)

# Normalised chemical species (:ce)
@parameters k
@variables E_1(t) F_1(t) q_1(t)
ce_dict = Dict(
    :description => "Chemical species (normalised)",
    :numports => 1,
    :params => Dict(
        k => "Biochemical Constant; exp(mu_0/RT)/V"
    ),
    :equations => [
        0 ~ log(k*q_1) - E_1,
        D(q_1) ~ F_1
    ]
)

# Chemical reaction (:Re)
@parameters r R T
@variables E_1(t) F_1(t) E_2(t) F_2(t)
Re_dict = Dict(
    :description => "Biochemical reaction",
    :numports => 2,
    :params => Dict(
        r => "Reaction Rate",
        R => "Universal Gas Constant",
        T => "Temperature"
    ),
    :equations => [
        0 ~ F_1 + F_2,
        0 ~ F_1 - r*(exp(E_1/R/T) - exp(E_2/R/T))
    ]
)

# Normalised chemical reaction (:re)
@parameters r R T
@variables E_1(t) F_1(t) E_2(t) F_2(t)
re_dict = Dict(
    :description => "Biochemical reaction (normalised)",
    :numports => 2,
    :params => Dict(
        r => "Reaction Rate"
    ),
    :equations => [
        0 ~ F_1 + F_2,
        0 ~ F_1 - r*(exp(E_1) - exp(E_2))
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