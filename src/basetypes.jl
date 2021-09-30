abstract type AbstractNode end

struct Component{N} <: AbstractNode
    type::Symbol
    name::AbstractString
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    function Component{N}(t::Symbol, n::AbstractString, np::Int, v::Int) where N
        new(t, n, ones(MVector{np,Bool}), Ref(v))
    end
end
Component(type::Symbol, name::AbstractString=string(type); numports::Int=1, vertex::Int=0) = 
    Component{numports}(type, name, numports, vertex)

struct Junction <: AbstractNode
    type::Symbol
    name::AbstractString
    vertex::RefValue{Int}
    Junction(t::Symbol, n::AbstractString=string(t); v::Int=0) = new(t, n, Ref(v))
end

struct Port 
    node::AbstractNode
    index::Int
    function Port(node::AbstractNode, index)
        ports = freeports(node)
        any(ports) || error("Node $node has no free ports")
        ports[index] || error("Port $index in node $node is already connected")
        new(node, index)
    end
end
Port(node::AbstractNode) = Port(node, nextfreeport(node))

struct Bond <: lg.AbstractSimpleEdge{Int}
    src::Port
    dst::Port
end
function Bond(srcnode::AbstractNode, dstnode::AbstractNode)
    Bond(Port(srcnode), Port(dstnode))
end

struct BondGraph <: lg.AbstractGraph{Int64}
    name::AbstractString
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(name::AbstractString="BG") = BondGraph(name, AbstractNode[], Bond[])

struct BondGraphNode <: AbstractNode
    bondgraph::BondGraph
    type::Symbol
    name::AbstractString
    freeports::Vector{Bool}
    vertex::RefValue{Int}
    function BondGraphNode(bg::BondGraph, type::Symbol=:BG, name::AbstractString=bg.name; vertex::Int=0)
        new(bg, type, name, Bool[], Ref(vertex))
    end
end