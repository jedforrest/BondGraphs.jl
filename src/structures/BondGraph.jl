# Type definition alias for MetaGraphs stored within Bond Graphs
# See https://github.com/JuliaGraphs/MetaGraphsNext.jl/issues/95#issuecomment-3572905564
const BondGraphMeta = MetaGraph{
    Int64,  # Code
    SimpleDiGraph{Int64},  # Graph
    Symbol,  # VertexLabel
    AbstractElement,  # VertexData
    Bond,  # EdgeData
    String  # GraphData
}

############################################################################################

"""
    BondGraph(name="BG") <: Graphs.AbstractGraph{Int64}
    BondGraph(name, nodes::Vector{AbstractNode}, bonds::Vector{Bond})

The bond graph object which contains a vector of `nodes` and `bonds`. All operations on
components or bonds must happen within the same bond graph. This inherits the methods of the
AbstractGraph type and so will work with

See also [`BondGraphNode`](@ref).
"""
mutable struct BondGraph <: AbstractGraph{Int}
    name::Symbol
    graph::MetaGraph
    sys::Union{System, Nothing}  # compiled later
    autocompile::Bool  # whether to recompile the System object whenever "system" is called
    function BondGraph(
            components::Vector{<:AbstractElement} = AbstractElement[],
            bonds::Vector{Bond} = Bond[];
            name = :BondGraph,
            autocompile = true
        )
        graph = metagraph(components, bonds, name)
        new(Symbol(name), graph, nothing, autocompile)
    end
end

# function BondGraph(name, bonds::Vector{Bond})
#     # TODO constuct from bonds only
#     BondGraph(Symbol(name), componentnames, bonds)
# end

function Base.show(io::IO, bg::BondGraph)
    print_str = "$(name(bg))"
    print_str *= isempty(bonds(bg)) ? "" : "\n  $(join(bonds(bg),"\n  "))"
    print(io, print_str)
end

############################################################################################

function metagraph(components::Vector{<:AbstractElement}, bonds::Vector{Bond}, bgname)
    graph = MetaGraph(
        DiGraph();
        label_type = Symbol,
        vertex_data_type = AbstractElement,
        edge_data_type = Bond,
        graph_data = string(bgname)
    )
    # check for duplicate names in components
    duplicate_names = [k for (k, v) in counter(name.(components)) if v > 1]
    if !isempty(duplicate_names)
        @warn "duplicate component names found: $(duplicate_names...)"
    end
    # vertices
    for comp in components
        graph[name(comp)] = comp
    end
    # edges
    for bond in bonds
        srcname, dstname = componentnames(bond)
        graph[srcname, dstname] = bond
    end
    graph
end

# TODO metagraph 'weight' can be the effort/flow pair in stored in the Bond

############################################################################################
name(bg::BondGraph) = bg.name

# components and edges are derived from the graph structure
components(bg::BondGraph) = [getindex(bg.graph, c) for c in labels(bg.graph)]
bonds(bg::BondGraph) = [getindex(bg.graph, s, d) for (s, d) in edge_labels(bg.graph)]

filterbytype(T::Type{<:AbstractElement}, vec) = filter(x -> x isa T, vec)
elements(bg::BondGraph) = filterbytype(BondElement, components(bg))
junctions(bg::BondGraph) = filterbytype(JunctionStructure, components(bg))

############################################################################################
# Extra simplification rules
log_exp_rules = [
    @rule(~a * log(~x) => log((~x) ^ (~a))),
    @acrule(~!a * log(~x) + ~!a * log(~y) => ~a * log(~x * ~y)),
    @acrule(~!a * exp(~!b * log(~x) + ~!c) => (~a) * ((~x)^(~b) + ~c)),
]
const rewriter = Postwalk(RestartedChain(log_exp_rules))

function system(bg::BondGraph; simplify = true)
    if isnothing(bg.sys) || bg.autocompile
        compile_system!(bg; simplify)
    end
    bg.sys
end

# MTK System compiler
function compile_system!(bg::BondGraph; simplify = true)
    comps = components(bg)
    subsyss = system.(comps)
    conn_eqns = connection_equation.(bonds(bg))

    basesys = System(conn_eqns, t, name = name(bg))
    sys = compose(basesys, subsyss...)
    if simplify && length(comps) > 0
        sys = structural_simplify(sys)
    end
    bg.sys = sys
end

function constitutive_relations(bg::BondGraph; sub_defaults = false)
    sys = system(bg)
    eqs = full_equations(sys)
    # for chemical exp/log cancelling
    # this only works with this syntax for some reason
    eqs = simplify.(simplify.(eqs); rewriter)
    if sub_defaults
        default_dict = defaults(sys)
        sub_dict = Dict(
            p => default_dict[p]
            for p in parameters(sys)
            if default_dict[p] isa Number
        )
        eqs = [substitute(eq, sub_dict) for eq in eqs]
    end
    eqs
end

# TODO display equations with variables indexed by the (graph) vertex of their component

############################################################################################
# Graph functions
# Most graph functions are passed on to the graph field within the bond graph struct

# eltype
eltype(::Type{BondGraph}) = Int
eltype(::BondGraph) = Int

# edgetype TODO check if correct
edgetype(::Type{BondGraph}) = SimpleEdge{Int}
edgetype(::BondGraph) = SimpleEdge{Int}

# edges
edges(bg::BondGraph) = edges(bg.graph)
has_edge(bg::BondGraph, src, dst) = has_edge(bg.graph, src, dst)
ne(bg::BondGraph) = ne(bg.graph)

# vertices
vertices(bg::BondGraph) = vertices(bg.graph)
has_vertex(bg::BondGraph, v) = has_vertex(bg.graph, v)
nv(bg::BondGraph) = nv(bg.graph)

# neighbors
inneighbors(bg::BondGraph, v) = inneighbors(bg.graph, v)
outneighbors(bg::BondGraph, v) = outneighbors(bg.graph, v)

# directed
is_directed(bg::Type{BondGraph}) = true
is_directed(bg::BondGraph) = true

# MetaGraph indexing
getindex(bg::BondGraph) = getindex(bg.graph)
getindex(bg::BondGraph, v) = getindex(bg.graph, v)
getindex(bg::BondGraph, s, d) = getindex(bg.graph, s, d)

setindex!(bg::BondGraph, data) = setindex!(bg.graph, data)
setindex!(bg::BondGraph, data, v) = setindex!(bg.graph, data, v)
setindex!(bg::BondGraph, data, s, d) = setindex!(bg.graph, data, s, d)

# graph mutations TODO

# add_vertex!
# rem_vertex!
# add_edge!
# rem_edge!

############################################################################################

inneighbor_comps(bg::BondGraph, elem::AbstractElement) = collect(inneighbor_labels(bg.graph, name(elem)))
outneighbor_comps(bg::BondGraph, elem::AbstractElement) = collect(outneighbor_labels(bg.graph, name(elem)))

function Base.in(elem::AbstractElement, bg::BondGraph)
    haskey(bg.graph, name(elem))
end

function Base.in(b::Bond, bg::BondGraph)
    srcname, dstname = componentnames(b)
    haskey(bg.graph, srcname, dstname)
end

############################################################################################

# # Easier referencing systems using a.b notation
# function getproperty(bg::BondGraph, sym::Symbol)
#     # Calling getfield explicitly avoids using "a.b" and causing a StackOverflowError
#     allnodes = getfield(bg, :nodes)
#     names = [getfield(n, :name) for n in allnodes]
#     symnodes = allnodes[names .== string(sym)]
#     if isempty(symnodes)
#         return getfield(bg, sym)
#     elseif length(symnodes) == 1
#         return symnodes[1]
#     else
#         return symnodes
#     end
# end
