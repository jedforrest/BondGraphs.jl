module BiochemicalLibrary

using ModelingToolkit

export biochemical_library

const _R = 8.314
const _T = 310.0

@parameters t
D = Differential(t)

@variables e[1:2](t) f[1:2](t)

# Chemical species (:Ce)
@parameters K R T
@variables q(t)
Ce_dict = Dict(
    :description => """
    Chemical species
    e = R*T*log(K*q)
    dq/dt = f
    K: Biochemical Constant; exp(mu_0/RT)/V [1.0],
    R: Universal Gas Constant [8.314],
    T: Temperature [310]
    q: Molar quantity [0.0]
    """,
    :numports => 1,
    :parameters => Dict(
        K => 1.0,
        R => _R,
        T => _T
    ),
    :states => Dict(
        q => 0.0
    ),
    :equations => [
        0 ~ R * T * log(K * q) - e[1],
        D(q) ~ f[1]
    ],
)

# Normalised chemical species (:ce)
@parameters K
@variables q(t)
ce_dict = Dict(
    :description => """
    Chemical species (normalised)
    e = log(K*q)
    dq/dt = f
    K: Biochemical Constant; exp(mu_0/RT)/V [1.0],
    q: Molar quantity [0.0]
    """,
    :numports => 1,
    :parameters => Dict(
        K => 1.0
    ),
    :states => Dict(
        q => 0.0
    ),
    :equations => [
        0 ~ log(K * q) - e[1],
        D(q) ~ f[1]
    ],
)

# Chemical reaction (:Re)
@parameters r R T
Re_dict = Dict(
    :description => """
    Biochemical reaction
    f₁ + f₂ = 0
    f₁ = r * [exp(e₁/RT) - exp(e₂/RT)]
    r: Reaction rate [1.0]
    R: Universal Gas Constant [8.314],
    T: Temperature [310]
    """,
    :numports => 2,
    :parameters => Dict(
        r => 1.0,
        R => _R,
        T => _T
    ),
    :equations => [
        0 ~ f[1] + f[2],
        0 ~ f[1] - r * (exp(e[1] / R / T) - exp(e[2] / R / T))
    ],
)

# Normalised chemical reaction (:re)
@parameters r
re_dict = Dict(
    :description => """
    Biochemical reaction (normalised)
    f₁ + f₂ = 0
    f₁ = r * [exp(e₁/RT) - exp(e₂/RT)]
    r: Reaction rate [1.0]
    """,
    :numports => 2,
    :parameters => Dict(
        r => 1.
    ),
    :equations => [
        0 ~ f[1] + f[2],
        0 ~ f[1] - r * (exp(e[1]) - exp(e[2]))
    ],
)

const biochemical_library = Dict(
    :Ce => Ce_dict,
    :ce => ce_dict,
    :Re => Re_dict,
    :re => re_dict
)

end