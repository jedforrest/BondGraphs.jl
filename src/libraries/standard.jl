module StandardLibrary

using ModelingToolkit

export standard_library

@parameters t
D = Differential(t)

@variables e[1:2](t) f[1:2](t)

# Linear resistance (:R)
@parameters R
R_dict = Dict(
    :description => """
    Generalised Linear Resistor
    e = R*f
    R: Resistance [1.0]
    """,
    :numports => 1,
    :parameters => Dict(
        R => 1.
    ),
    :equations => [0 ~ e[1] - R * f[1]]
)

# Linear capacitor (:C)
@parameters C
@variables q(t)
C_dict = Dict(
    :description => """
    Generalised Linear Capacitor
    e = (1/C)*q
    dq/dt = f
    C: Capacitance [1.0]
    q: Generalised position [0.0]
    """,
    :numports => 1,
    :parameters => Dict(
        C => 1.
    ),
    :states => Dict(
        q => 0.
    ),
    :equations => [
        0 ~ q / C - e[1],
        D(q) ~ f[1]
    ],
)

# Linear inductance (:I)
@parameters L
@variables p(t)
I_dict = Dict(
    :description => """
    Generalised Linear Inductor
    f = (1/L)*p
    dp/dt = e
    L: Inductance [1.0]
    p: Generalised momentum [0.0]
    """,
    :numports => 1,
    :parameters => Dict(
        L => 1.
    ),
    :states => Dict(
        p => 0.
    ),
    :equations => [
        0 ~ p / L - f[1],
        D(p) ~ e[1]
    ],
)

# Source of effort (:Se)
@parameters E
Se_dict = Dict(
    :description => """
    Effort Source
    e = E
    E: Effort (source) [1.0]
    """,
    :numports => 1,
    :parameters => Dict(
        E => 1.
    ),
    :equations => [0 ~ E - e[1]],
)

# Source of flow (:Sf)
@parameters F
Sf_dict = Dict(
    :description => """
    Flow Source
    f = F
    F: Flow (source) [1.0]
    """,
    :numports => 1,
    :parameters => Dict(
        F => 1.
    ),
    :equations => [0 ~ F - f[1]],
)

# Transformer (:TF)
@parameters n
TF_dict = Dict(
    :description => """
    Linear Transformer
    e₂ = n*e₁
    f₁ = n*f₂
    n = Winding ratio [1.0]
    """,
    :numports => 2,
    :parameters => Dict(
        n => 1.
    ),
    :equations => [
        0 ~ e[2] - n * e[1],
        0 ~ f[1] - n * f[2]
    ],
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