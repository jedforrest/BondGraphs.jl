# Elements define a particular type of bond graph component (e.g. an electrical capacitor).
# An element will have constitutive relations (equations), efforts/flows, and any states.
# A bond graph can have multiple instances of an element.
# See Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

############################################################
abstract type AbstractElement end
abstract type BondElement <: AbstractElement end

"""`R` component"""
@kwdef struct Resistance <: BondElement
    symbol::Symbol = :R
end

#########################
abstract type StorageElement <: BondElement end

"""`C` component"""
@kwdef struct Compliance <: StorageElement
    symbol::Symbol = :C
end

"""`I` component"""
@kwdef struct Inertance <: StorageElement
    symbol::Symbol = :I
end

############################################################
abstract type SourceElement <: BondElement end

"""`Se` component"""
@kwdef struct EffortSource <: SourceElement
    symbol::Symbol = :Se
end

"""`Sf` component"""
@kwdef struct FlowSource <: SourceElement
    symbol::Symbol = :Sf
end

"""`SS` component"""
@kwdef struct SourceSensor <: SourceElement
    symbol::Symbol = :SS
end

############################################################
abstract type JunctionStructure <: AbstractElement end
abstract type ParametricJunction <: JunctionStructure end

# Junctions are a power-conserving transformation between the effort/flows of each port
# They have defined relationships between the inputs and outputs
# New junction types would require a new struct definition

"""`TF` component"""
@kwdef struct Transformer <: ParametricJunction
    symbol::Symbol = :TF
end

"""`GY` component"""
@kwdef struct Gyrator <: ParametricJunction
    symbol::Symbol = :GY
end

############################################################
abstract type NonParametricJunction <: JunctionStructure end

# Since the number of ports is unknown and can change,
# we instead store functions that map efforts and flows to equations
"""`0`-junction"""
@kwdef struct ZeroJunction <: NonParametricJunction
    symbol::Symbol = :𝟎
end

"""`1`-junction"""
@kwdef struct OneJunction <: NonParametricJunction
    symbol::Symbol = :𝟏
end

# compatibility
const EqualEffort = ZeroJunction
const EqualFlow = OneJunction

############################################################

symbol(elem::AbstractElement) = elem.symbol

function Base.show(io::IO, elem::AbstractElement)
    print(io, "$(symbol(elem))::$(typeof(elem))")
end
