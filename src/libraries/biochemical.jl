module BiochemicalLibrary

using ModelingToolkit
using OrderedCollections

export biochemical_library

@parameters t
D = Differential(t)

# Biochemical constant, reaction rate, Universal gas constant, Temperature
@parameters k, r, R, T

# Molar Quantity
@variables q(t)

@variables E[1:2](t) F[1:2](t)

# Chemical species (:Ce)
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
        0 ~ R * T * log(k * q) - E[1],
        D(q) ~ F[1]
    ]
)

# Normalised chemical species (:ce)
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
        0 ~ log(k * q) - E[1],
        D(q) ~ F[1]
    ]
)

# Chemical reaction (:Re)
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
        0 ~ F[1] - r * (exp(E[1] / R / T) - exp(E[2] / R / T))
    ]
)

# Normalised chemical reaction (:re)
re_dict = Dict(
    :description => "Biochemical reaction (normalised)",
    :numports => 2,
    :parameters => OrderedDict(
        r => "Reaction Rate"
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - r * (exp(E[1]) - exp(E[2]))
    ]
)

const biochemical_library = Dict(
    :Ce => Ce_dict,
    :ce => ce_dict,
    :Re => Re_dict,
    :re => re_dict
)

end