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
Component(metamodel::Symbol; name::String=string(metamodel), numports::Int=1, vertex::Int=0) = 
    Component{numports}(metamodel, name, numports, vertex)

struct Junction <: AbstractNode
    metamodel::Symbol
    vertex::RefValue{Int}
    Junction(m::Symbol; v::Int=0) = new(m, Ref(v))
end

struct Port 
    node::AbstractNode
    index::Int
    function Port(node::AbstractNode, index=nothing)
        ports = freeports(node)
        any(ports) || error("Node $node has no free ports")
        if isnothing(index)
            index = findfirst(ports)
        else
            ports[index] || error("Port $index in node $node is already connected")
        end
        new(node, index)
    end
end

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
BondGraph(metamodel::Symbol=:BG; name::String="BG") = BondGraph(metamodel, name, AbstractNode[], Bond[])

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Ports
freeports(n::Component) = n.freeports
freeports(n::Junction) = [true]
numports(n::Component) = length(n.freeports)
numports(n::Junction) = Inf
updateport!(n::AbstractNode, idx::Int) = freeports(n)[idx] = !freeports(n)[idx]

# Nodes in Bon`ds
srcnode(b::Bond) = b.src.node
dstnode(b::Bond) = b.dst.node
in(n::AbstractNode, b::Bond) = n == srcnode(b) || n == dstnode(b)

# I/O
show(io::IO, node::Component) = print(io, "$(node.metamodel):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.metamodel)")
show(io::IO, b::Bond) = print(io, "Bond $(srcnode(b)) ⇀ $(dstnode(b))")
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.metamodel):$(bg.name) ($(lg.nv(bg)) Nodes, $(lg.ne(bg)) Bonds)")
