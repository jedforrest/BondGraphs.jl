# Elements define a particular type of bond graph component (e.g. an electrical capacitor).
# An element will have constitutive relations (equations), efforts/flows, and any states.
# A bond graph can have multiple instances of an element.
# See Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

############################################################
abstract type AbstractElement end

abstract type BondElement <: AbstractElement end
abstract type StorageElement <: BondElement end
abstract type SourceElement <: BondElement end

"""`C` component"""
struct StaticStorageElement <: StorageElement
    sys::System
    efforts::Vector
    flows::Vector
    states::Vector
    function StaticStorageElement(sys::System, efforts, flows, states)
        check_num_ports(efforts, flows)
        new(sys, efforts, flows, states)
    end
end

"""`I` component"""
struct DynamicStorageElement <: StorageElement
    sys::System
    efforts::Vector
    flows::Vector
    states::Vector
    function DynamicStorageElement(sys::System, efforts, flows, states)
        check_num_ports(efforts, flows)
        new(sys, efforts, flows, states)
    end
end

"""`R` component"""
struct DissipatorElement <: BondElement
    sys::System
    efforts::Vector
    flows::Vector
    function DissipatorElement(sys::System, efforts, flows)
        check_num_ports(efforts, flows)
        new(sys, efforts, flows)
    end
end

############################################################
"""`Se` component"""
struct EffortSource <: SourceElement
    sys::System
    efforts::Vector
    flows::Vector
    function EffortSource(sys::System, efforts, flows)
        numports = length(efforts)
        numports == 1 || error("Must have exactly 1 port ($numports)")
        new(sys, efforts, flows)
    end
end

"""`Sf` component"""
struct FlowSource <: SourceElement
    sys::System
    efforts::Vector
    flows::Vector
    function FlowSource(sys::System, efforts, flows)
        numports = length(flows)
        numports == 1 || error("Must have exactly 1 port ($numports)")
        new(sys, efforts, flows)
    end
end

"""`SS` component"""
struct SourceSensor <: SourceElement
    sys::System
    efforts::Vector
    flows::Vector
    function SourceSensor(sys::System, efforts, flows)
        numports = length(efforts)
        (numports == length(flows) == 1) || error("Must have exactly 1 port ($numports)")
        new(sys, efforts, flows)
    end
end

############################################################
abstract type JunctionStructure <: AbstractElement end
abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

"""`TF` component"""
struct Transformer <: ParametricJunction
    sys::System
    efforts::Vector
    flows::Vector
    function Transformer(sys::System, efforts, flows)
        numports = length(efforts)
        numports >= 2 ||
            error("Transformer must have at least 2 ports ($numports ports given)")
        new(sys, efforts, flows)
    end
end
"""`GY` component"""
struct Gyrator <: ParametricJunction
    sys::System
    efforts::Vector
    flows::Vector
    function Gyrator(sys::System, efforts, flows)
        numports = length(efforts)
        numports >= 2 || error("Gyrator must have at least 2 ports ($numports ports given)")
        new(sys, efforts, flows)
    end
end

############################################################
# Since the number of ports is unknown and can change,
# we instead store functions that map efforts/flows to equations
"""`0`-junction"""
struct EqualEffort <: NonParametricJunction
    effort_fn::Function
    flow_fn::Function
    function EqualEffort()
        effort_fn(es) = [es[1] ~ e for e in es[2:end]]
        flow_fn(fs) = sum(fs) ~ 0
        new(effort_fn, flow_fn)
    end
end
"""`1`-junction"""
struct EqualFlow <: NonParametricJunction
    effort_fn::Function
    flow_fn::Function
    function EqualFlow()
        effort_fn(es) = sum(es) ~ 0
        flow_fn(fs) = [fs[1] ~ f for f in fs[2:end]]
        new(effort_fn, flow_fn)
    end
end
############################################################

""" Generic BondGraph Element constructor from a vector of equations"""
function (E::Type{<:AbstractElement})(eqs::Vector{Equation}, efforts, flows, states=[])
    sys = System(eqs, t; name=Symbol(E))
    if hasfield(E, :states)
        return E(sys, efforts, flows, states)
    else
        return E(sys, efforts, flows)
    end
end

function check_num_ports(efforts, flows)
    numports = length(efforts)
    numports >= 1 || error("Must have at least 1 port ($numports)")
    (numports == length(flows)) ||
        error("Number of efforts and flows don't match ($numports)")
    return nothing
end

equations(elem::AbstractElement) = equations(elem.sys)
equations(::NonParametricJunction) = nothing  # FIXME

numports(elem::AbstractElement) = length(elem.efforts)
numports(::SourceElement) = 1
numports(::NonParametricJunction) = Inf

# Used when displaying in a graph.
# TODO these can just be included in the struct definitions above (kwdef)
glyph(::DissipatorElement) = :R
glyph(::StaticStorageElement) = :C
glyph(::DynamicStorageElement) = :I
glyph(::EffortSource) = :Se
glyph(::FlowSource) = :Sf
glyph(::Transformer) = :TF
glyph(::Gyrator) = :GY
glyph(::JunctionStructure) = :J
glyph(::EqualEffort) = :𝟎
glyph(::EqualFlow) = :𝟏

function Base.show(io::IO, elem::T) where {T <: AbstractElement}
    eqs = equations(elem)
    print_str = "$T{$(numports(elem))}"
    print_str *= isempty(eqs) ? "" : "\n  $(join(eqs,"\n  "))"
    print(io, print_str)
end
Base.show(io::IO, ::T) where {T <: NonParametricJunction} = print(io, T)
