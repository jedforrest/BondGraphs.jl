# Elements define a particular type of bond graph component (e.g. an electrical capacitor).
# An element will have constitutive relations (equations), efforts/flows, and any states.
# A bond graph can have multiple instances of an element.
# See Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

############################################################
abstract type AbstractElement end
abstract type BondElement <: AbstractElement end

"""`R` component"""
@kwdef struct DissipatorElement <: BondElement
    symbol::Symbol = :R
    numports::Int = 1
end

#########################
abstract type StorageElement <: BondElement end

"""`C` component"""
@kwdef struct StaticStorageElement <: StorageElement
    symbol::Symbol = :C
    numports::Int = 1
end

"""`I` component"""
@kwdef struct DynamicStorageElement <: StorageElement
    symbol::Symbol = :I
    numports::Int = 1
end

############################################################
abstract type SourceElement <: BondElement end

"""`Se` component"""
@kwdef struct EffortSource <: SourceElement
    symbol::Symbol = :Se
    numports::Int = 1
end

"""`Sf` component"""
@kwdef struct FlowSource <: SourceElement
    symbol::Symbol = :Sf
    numports::Int = 1
end

"""`SS` component"""
@kwdef struct SourceSensor <: SourceElement
    symbol::Symbol = :SS
    numports::Int = 1
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
    numports::Int = 2
    # function Transformer(n::Real; name)
    #     ports = [Port(:port_1, name), Port(:port_2, name)]
    #     e1, e2 = effort.(ports)
    #     f1, f2 = flow.(ports)
    #     eqs = [
    #         e2 ~ n * e1
    #         f1 ~ -n * f2
    #     ]
    #     sys = System(eqs, t; name, systems=[ports[1].sys, ports[2].sys])
    #     new(name, sys, n, ports)
    # end
end

"""`GY` component"""
@kwdef struct Gyrator <: ParametricJunction
    symbol::Symbol = :GY
    numports::Int = 2
    # function Gyrator(r::Real; name)
    #     ports = [Port(:port_1, name), Port(:port_2, name)]
    #     e1, e2 = effort.(ports)
    #     f1, f2 = flow.(ports)
    #     eqs = [
    #         e2 ~ r * f1
    #         e1 ~ -r * f2
    #     ]
    #     sys = System(eqs, t; name, systems=[ports[1].sys, ports[2].sys])
    #     new(name, sys, r, ports)
    # end
end

############################################################
abstract type NonParametricJunction <: JunctionStructure end

# Since the number of ports is unknown and can change,
# we instead store functions that map efforts and flows to equations
"""`0`-junction"""
@kwdef struct ZeroJunction <: NonParametricJunction
    symbol::Symbol = :𝟎
    # function EqualEffort(; name)
    #     sys = System(Equation[], t; name)
    #     new(name, sys, [Port(:port_1, name)])
    # end
end

"""`1`-junction"""
@kwdef struct OneJunction <: NonParametricJunction
    symbol::Symbol = :𝟏
    # function EqualFlow(; name)
    #     sys = System(Equation[], t; name)
    #     new(name, sys, [Port(:port_1, name)])
    # end
end

# compatibility
const EqualEffort = ZeroJunction
const EqualFlow = OneJunction

############################################################

symbol(elem::AbstractElement) = elem.symbol

numports(elem::BondElement) = elem.numports
numports(junc::ParametricJunction) = junc.numports
numports(::NonParametricJunction) = Inf

############################################################
# Overloading Base

function Base.show(io::IO, elem::AbstractElement)
    print(io, "$(typeof(elem))($(symbol(elem)), $(numports(elem)))")
end

function Base.show(io::IO, elem::NonParametricJunction)
    print(io, "$(typeof(elem))($(symbol(elem)))")
end
