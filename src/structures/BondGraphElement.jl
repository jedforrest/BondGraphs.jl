"""
    BondGraphElement(bg::BondGraph, name=name(bg); deep_copy=false)

Convert a `BondGraph` into a component that can be added in another level bond graph.
Componets can be exposed to the outer bond graph by replacing them with a [`SourceSensor`](@ref)
type using the [`swap!`](@ref) function.

See also [`BondGraph`](@ref).
"""
struct BondGraphElement <: AbstractElement
    name::Symbol
    bondgraph::BondGraph
    # sys::System
    efforts::Vector
    flows::Vector
    ports::Vector{Port}
end
function BondGraphElement(bg::BondGraph; name = name(bg), deep_copy = false)
    _bg = deep_copy ? deepcopy(bg) : bg

    sources = filterbytype(SourceSensor, components(_bg))
    isempty(sources) && @warn("$bg has no exposed ports (SourceSensors)")

    es = efforts.(sources)
    fs = flows.(sources)
    ports = Port[Port(ss.name, name) for ss in sources]

    BondGraphElement(name, _bg, es, fs, ports)
end

# Easier referencing systems using a.b notation
# function getproperty(bgn::BondGraphElement, sym::Symbol)
#     bg = getfield(bgn, :bondgraph)
#     try
#         return getproperty(bg, sym)
#     catch
#         return getfield(bgn, sym)
#     end
# end

icon(::BondGraphElement) = :BG

function system(bge::BondGraphElement)

    sys = system(bge.bondgraph)

end
