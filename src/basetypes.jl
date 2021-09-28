abstract type AbstractNode end

struct Component{N} <: AbstractNode
    metamodel::Symbol
    name::AbstractString
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    function Component{N}(m::Symbol, n::AbstractString, np::Int, v::Int) where N
        new(m, n, ones(MVector{np,Bool}), Ref(v))
    end
end
Component(metamodel::Symbol, name::AbstractString=string(metamodel); numports::Int=1, vertex::Int=0) = 
    Component{numports}(metamodel, name, numports, vertex)

struct Junction <: AbstractNode
    metamodel::Symbol
    name::AbstractString
    vertex::RefValue{Int}
    Junction(m::Symbol, n::AbstractString=string(m); v::Int=0) = new(m, n, Ref(v))
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
    metamodel::Symbol
    name::AbstractString
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(metamodel::Symbol=:BG, name::AbstractString=string(metamodel)) = BondGraph(metamodel, name, AbstractNode[], Bond[])
BondGraph(name::AbstractString) = BondGraph(:BG, name)
