"""
    BondGraph(name="BG") <: Graphs.AbstractGraph{Int64}
    BondGraph(name, nodes::Vector{AbstractNode}, bonds::Vector{Bond})

The bond graph object which contains a vector of `nodes` and `bonds`. All operations on
components or bonds must happen within the same bond graph. This inherits the methods of the
AbstractGraph type and so will work with

See also [`BondGraphNode`](@ref).
"""
struct BondGraph <: g.AbstractGraph{Int64}
    name::AbstractString
    nodes::Vector{T} where {T<:AbstractNode}
    bonds::Vector{Bond}
end
function BondGraph(name="BG")
    BondGraph(string(name), AbstractNode[], Bond[])
end

# PROPERTIES
name(bg::BondGraph) = bg.name

nodes(bg::BondGraph) = bg.nodes

bonds(bg::BondGraph) = bg.bonds

# AbstractNode properties
function _nested_bg_variables(bg::BondGraph, var_function::Function)
    OrderedDict(comp => Dict(var for var in var_function(comp)) for comp in components(bg))
end

parameters(bg::BondGraph) = _nested_bg_variables(bg, parameters)

globals(bg::BondGraph) = _nested_bg_variables(bg, globals)

states(bg::BondGraph) = _nested_bg_variables(bg, states)

controls(bg::BondGraph) = _nested_bg_variables(bg, controls)

all_variables(bg::BondGraph) = _nested_bg_variables(bg, all_variables)

function equations(bg::BondGraph; simplify_eqs=true)
    isempty(bg.nodes) && return Equation[]
    sys = ODESystem(bg; simplify_eqs)
    return equations(sys)
end

has_controls(bg::BondGraph) = any(.!isempty.(controls.(nodes(bg))))

# Filtering
components(bg::BondGraph) = filter(x -> x isa Component, bg.nodes)
junctions(bg::BondGraph) = filter(x -> x isa Junction, bg.nodes)

"""
    getnodes(bg::BondGraph, type)

Return all nodes a particular bond graph type in the bond graph `bg`.

`type` can be a DataType (e.g. Component{1}), a string (e.g. "C"), or a vector of strings.
"""
getnodes(bg::BondGraph, T::DataType) = filter(n -> n isa T, bg.nodes)
getnodes(bg::BondGraph, t::AbstractString) = filter(n -> "$(type(n)):$(name(n))" == t, bg.nodes)
getnodes(bg::BondGraph, ts::Vector{T} where T <: AbstractString) = vcat((getnodes(bg, t) for t in ts)...)

"""
    getbonds(bg::BondGraph, n1::AbstractNode, n2::AbstractNode)
    getbonds(bg::BondGraph, (n1, n2))

Return the bond in `bg` connecting nodes `n1` and `n2`, if it exists.
"""
getbonds(bg::BondGraph, t::Tuple) = getbonds(bg, t[1], t[2])
getbonds(bg::BondGraph, n1::AbstractNode, n2::AbstractNode) = filter(b -> n1 in b && n2 in b, bg.bonds)


# Base functions
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.name) ($(g.nv(bg)) Nodes, $(g.ne(bg)) Bonds)")

# Easier referencing systems using a.b notation
function getproperty(bg::BondGraph, sym::Symbol)
    # Calling getfield explicitly avoids using "a.b" and causing a StackOverflowError
    allnodes = getfield(bg, :nodes)
    names = [getfield(n, :name) for n in allnodes]
    symnodes = allnodes[names.==string(sym)]
    if isempty(symnodes)
        return getfield(bg, sym)
    elseif length(symnodes) == 1
        return symnodes[1]
    else
        return symnodes
    end
end

# Conversion to common graph types
g.SimpleGraph(bg::BondGraph) = g.SimpleGraph(g.SimpleDiGraph(bg))
g.SimpleDiGraph(bg::BondGraph) = g.SimpleDiGraph(g.adjacency_matrix(bg))

"""
    BondGraphNode(bg::BondGraph, name=name(bg); deep_copy=false)

Convert a `BondGraph` into a component that can be added in another level bond graph.
Componets can be exposed to the outer bond graph by replacing them with a [`SourceSensor`](@ref)
type using the [`swap!`](@ref) function.

See also [`BondGraph`](@ref).
"""
struct BondGraphNode <: AbstractNode
    bondgraph::BondGraph
    type::AbstractString
    name::AbstractString
    ports::OrderedDict{Any,Bool}
    vertex::RefValue{Int}
end
function BondGraphNode(bg::BondGraph, name=name(bg); vertex::Int=0, deep_copy=false)
    _bg = deep_copy ? deepcopy(bg) : bg

    exposed_ports = getnodes(_bg, SourceSensor)
    ports = Dict(i => false for i in 1:length(exposed_ports))

    BondGraphNode(_bg, "BG", string(name), ports, Ref(vertex))
end

# Easier referencing systems using a.b notation
function getproperty(bgn::BondGraphNode, sym::Symbol)
    bg = getfield(bgn, :bondgraph)
    try
        return getproperty(bg, sym)
    catch
        return getfield(bgn, sym)
    end
end

exposed(bgn::BondGraphNode) = getnodes(bgn.bondgraph, SourceSensor)
