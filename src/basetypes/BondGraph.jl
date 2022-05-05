# BONDGRAPH
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
parameters(bg::BondGraph) = Dict(p for c in components(bg) for p in parameters(c))

globals(bg::BondGraph) = Dict(g for c in components(bg) for g in globals(c))

states(bg::BondGraph) = Dict(x for c in components(bg) for x in states(c))

controls(bg::BondGraph) = Dict(u for c in components(bg) for u in controls(c))

all_variables(bg::BondGraph) = Dict(v for c in components(bg) for v in all_variables(c))

function equations(bg::BondGraph; simplify_eqs=true)
    isempty(bg.nodes) && return Equation[]
    sys = ODESystem(bg; simplify_eqs=simplify_eqs)
    return equations(sys)
end

# Filtering
components(bg::BondGraph) = filter(x -> x isa Component, bg.nodes)
junctions(bg::BondGraph) = filter(x -> x isa Junction, bg.nodes)

getnodes(bg::BondGraph, T::DataType) = filter(x -> x isa T, bg.nodes)
getnodes(bg::BondGraph, t::AbstractString) = filter(x -> type(x) == t, bg.nodes)

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


# BONDGRAPH NODE
struct BondGraphNode <: AbstractNode
    bondgraph::BondGraph
    type::AbstractString
    name::AbstractString
    exposed::Vector{SourceSensor}
    freeports::Vector{Bool}
    vertex::RefValue{Int}
end
function BondGraphNode(bg::BondGraph, name=name(bg); vertex::Int=0, deepcopy=false)
    _bg = deepcopy ? deepcopy(bg) : bg

    exposed_ports = getnodes(_bg, SourceSensor)
    freeports = fill(true, length(exposed_ports))

    BondGraphNode(_bg, "BG", string(name), exposed_ports, freeports, Ref(vertex))
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

exposed(bgn::BondGraphNode) = bgn.exposed
