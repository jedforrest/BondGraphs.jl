# BONDGRAPH
struct BondGraph <: g.AbstractGraph{Int64}
    type::Symbol
    name::Symbol
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(name=:BG) = BondGraph(:BG, Symbol(name), AbstractNode[], Bond[])

# BONDGRAPH NODE
# struct BondGraphNode <: AbstractNode
#     bondgraph::BondGraph
#     type::Symbol
#     name::Symbol
#     freeports::Vector{Bool}
#     vertex::RefValue{Int}
#     function BondGraphNode(bg::BondGraph, type::Symbol=:BG, name::Symbol=bg.name; vertex::Int=0)
#         new(bg, type, name, Bool[], Ref(vertex))
#     end
# end



# PROPERTIES
type(bg::BondGraph) = bg.type

name(bg::BondGraph) = bg.name

nodes(bg::BondGraph) = bg.nodes

bonds(bg::BondGraph) = bg.bonds


# State variables
# function state_vars(m::BondGraph)
#     i = 0
#     states = Vector{Tuple{Union{AbstractNode,BondGraph},Num}}([])
#     for c in components(m), x in state_vars(c)
#         push!(states,(c,x))
#     end
#     @variables x[1:length(states)](t)
#     return OrderedDict(x[i] => y for (i,y) in enumerate(states))
# end


# Filtering
components(bg::BondGraph) = filter(x -> x isa Component, bg.nodes)
junctions(bg::BondGraph) = filter(x -> x isa Junction, bg.nodes)

getnodes(bg::BondGraph, t::Symbol) = filter(x -> type(x) == t, bg.nodes)
getnodes(bg::BondGraph, n) = filter(x -> name(x) == Symbol(n), bg.nodes)

getbonds(bg::BondGraph, t::Tuple) = getbonds(bg, t[1], t[2])
getbonds(bg::BondGraph, n1::AbstractNode, n2::AbstractNode) = filter(b -> n1 in b && n2 in b, bg.bonds)


# Base functions
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.type):$(bg.name) ($(g.nv(bg)) Nodes, $(g.ne(bg)) Bonds)")

# Easier referencing systems using a.b notation
# TODO create tests
function getproperty(bg::BondGraph, sym::Symbol)
    # Calling getfield explicitly avoids using "a.b" and causing a StackOverflowError
    allnodes = getfield(bg, :nodes)
    names = [getfield(n, :name) for n in allnodes]
    symnodes = allnodes[names .== sym]
    if isempty(symnodes)
        return getfield(bg, sym)
    elseif length(symnodes) == 1
        return symnodes[1]
    else
        return symnodes
    end
end

# function getproperty(bgn::BondGraphNode, sym::Symbol)
#     bg = getfield(bgn, :bondgraph)
#     try 
#         return getproperty(bg, sym)
#     catch
#         return getfield(bgn, sym)
#     end
# end