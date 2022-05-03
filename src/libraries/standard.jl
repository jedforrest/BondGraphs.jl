module StandardLibrary

using ModelingToolkit

export standard_library

# Schema Explanation
# description -> written description of the component: name, equation, parameter and variable definitions
# numports -> the number of ports the component has (fixed)
# variables ->
#   parameters -> constant parameters, unique for each component instance (will have a namespace)
#   globals -> global parameters (no namespace in the model)
#   states -> time-dependent state variables
#   controls -> time-dependent parameters, can accept an arbitrary julia function (can also remain constant)
# equations -> symbolic description of the constitutive equations

@parameters t
D = Differential(t)

@variables E[1:2](t) F[1:2](t)

# Linear resistance (:R)
@parameters R
R_dict = Dict(
    :description=>"""
    Generalised Linear Resistor
    e = R*f
    R: Resistance [1.0]
    """,
    :numports => 1,
    :variables => Dict(
        :parameters => Dict(
            R => 1.0
        ),
    ),
    :equations=>[0 ~ E[1] - R * F[1]]
)

# Linear capacitor (:C)
@parameters C
@variables q(t)
C_dict = Dict(
    :description=>"""
    Generalised Linear Capacitor
    e = (1/C)*q
    dq/dt = f
    C: Capacitance [1.0]
    q: Generalised position [0.0]
    """,
    :numports => 1,
    :variables => Dict(
        :parameters => Dict(
            C => 1.0
        ),
        :states => Dict(
            q => 0.0
        ),
    ),
    :equations => [
        0 ~ q / C - E[1],
        D(q) ~ F[1]
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
    :variables => Dict(
        :parameters => Dict(
            L => 1.0
        ),
        :states => Dict(
            p => 0.0
        ),
    ),
    :equations=>[
        0 ~ p / L - F[1],
        D(p) ~ E[1]
    ],
)

# Source of effort (:Se)
@parameters es(t)
Se_dict = Dict(
    :description => """
      Effort Source
      e = eₛ
      eₛ: Effort (source) [1.0]
      """,
    :numports => 1,
    :variables => Dict(
        :controls => Dict(
            es => (t -> 1.0)
        ),
    ),
    :equations=>[0 ~ es - E[1]],
)

# Source of flow (:Sf)
@parameters fs(t)
Sf_dict = Dict(
    :description => """
      Flow Source
      f = fₛ
      fₛ: Flow (source) [1.0]
      """,
    :numports => 1,
    :variables => Dict(
        :controls => Dict(
            fs => (t -> 1.0)
        ),
    ),
    :equations=>[0 ~ fs + F[1]],
)

# Transformer (:TF)
@parameters n
TF_dict = Dict(
    :description=>"""
    Linear Transformer
    e₂ = n*e₁
    f₁ = n*f₂
    n = Winding ratio [1.0]
    """,
    :numports => 2,
    :variables => Dict(
        :parameters => Dict(
            n => 1.0
        ),
    ),
    :equations=>[
        0 ~ E[2] - n * E[1],
        0 ~ F[1] - n * F[2]
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
