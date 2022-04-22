module BiochemicalLibrary

using ModelingToolkit

export biochemical_library

const _R = 8.314
const _T = 310.0

# GlobalScope parameters do not have a namespace
@parameters R T
R, T = GlobalScope(R), GlobalScope(T)

@parameters t
D = Differential(t)

@variables E[1:2](t) F[1:2](t)

# Chemical species (:Ce)
@parameters K #R T
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
        K => 1.0
    ),
    :globals => Dict(
        R => _R,
        T => _T
    ),
    :states => Dict(
        q => 0.0
    ),
    :equations => [
        0 ~ R * T * log(K * q) - E[1],
        D(q) ~ F[1]
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
        0 ~ log(K * q) - E[1],
        D(q) ~ F[1]
    ],
)

# Chemical reaction (:Re)
@parameters r #R T
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
        r => 1.0
    ),
    :globals => Dict(
        R => _R,
        T => _T
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - r * (exp(E[1] / R / T) - exp(E[2] / R / T))
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
        r => 1.0
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - r * (exp(E[1]) - exp(E[2]))
    ],
)

const biochemical_library = Dict(
    :Ce => Ce_dict,
    :ce => ce_dict,
    :Re => Re_dict,
    :re => re_dict
)

end