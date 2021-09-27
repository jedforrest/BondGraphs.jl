abstract type AbstractNode end

struct Component{N} <: AbstractNode
    type::Symbol
    name::AbstractString
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    parameters::Vector{Num}
    state_vars::Vector{Num}
    equations::Vector{Equation}
    function Component{N}(m::Symbol, n::AbstractString, np::Int, v::Int,
         p::Vector, x::Vector, eq::Vector) where N
        new(m, n, ones(MVector{np,Bool}), Ref(v), p, x, eq)
    end
end
function Component(type::Symbol, name::String=string(type);
    numports::Int=1, vertex::Int=0, parameters::Vector=[], state_vars::Vector=[],
    equations::Vector=[])
    Component{numports}(type, name, numports, vertex, parameters, state_vars, equations)
end

struct Junction <: AbstractNode
    type::Symbol
    name::AbstractString
    vertex::RefValue{Int}
    Junction(m::Symbol; n::AbstractString=string(m), v::Int=0) = new(m, n, Ref(v))
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
    type::Symbol
    name::AbstractString
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(type::Symbol=:BG, name::AbstractString=string(type)) = BondGraph(type, name, AbstractNode[], Bond[])
BondGraph(name::AbstractString) = BondGraph(:BG, name)

# New component
function new(type,name::String=string(type);library=standard_library)
    d = library[type]
    p = collect(keys(d[:parameters]))
    if haskey(d,:state_vars)
        x = collect(keys(d[:state_vars]))
    else
        x = Vector{Num}([])
    end
    Component(type, name; numports=d[:numports], parameters=p, state_vars=x, equations=d[:equations])
end

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Ports
freeports(n::Component) = n.freeports
freeports(n::Junction) = [true]
numports(n::Component) = length(n.freeports)
numports(n::Junction) = Inf
updateport!(n::AbstractNode, idx::Int) = freeports(n)[idx] = !freeports(n)[idx]
nextfreeport(n::AbstractNode) = findfirst(freeports(n))

# Nodes in Bonds
srcnode(b::Bond) = b.src.node
dstnode(b::Bond) = b.dst.node
in(n::AbstractNode, b::Bond) = n == srcnode(b) || n == dstnode(b)

# Equations
equations(n::Component) = n.equations

# Parameters
params(n::Component) = n.parameters

# State variables
state_vars(n::Component) = n.state_vars

# I/O
show(io::IO, node::Component) = print(io, "$(node.type):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.type)")
show(io::IO, port::Port) = print(io, "Port $(port.node) ($(port.index))")
show(io::IO, b::Bond) = print(io, "Bond $(srcnode(b)) ⇀ $(dstnode(b))")
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(bg.type):$(bg.name) ($(lg.nv(bg)) Nodes, $(lg.ne(bg)) Bonds)")
