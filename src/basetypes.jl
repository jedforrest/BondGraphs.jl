abstract type AbstractNode end

struct Component{N} <: AbstractNode
    type::Symbol
    name::Symbol
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    parameters::Vector{Num}
    state_vars::Vector{Num}
    default::Dict{Num,Any}
    equations::Vector{Equation}
    function Component{N}(m::Symbol, n::Symbol, np::Int, v::Int,
         p::Vector, x::Vector, d::Dict, eq::Vector) where N
        new(m, n, ones(MVector{np,Bool}), Ref(v), p, x, d, eq)
    end
end
function Component(type::Symbol, name::Symbol=type;
    numports::Int=1, vertex::Int=0, parameters::Vector=[], state_vars::Vector=[],
    default::Dict=Dict(), equations::Vector=[])
    Component{numports}(type, name, numports, vertex, parameters, state_vars, default, equations)
end

abstract type Junction <: AbstractNode end

struct EqualEffort <: Junction
    name::Symbol
    freeports::Vector{Bool}
    vertex::RefValue{Int}
    EqualEffort(; name::Symbol=Symbol("0"), v::Int=0) = new(name, [], Ref(v))
end

struct EqualFlow <: Junction
    name::Symbol
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    EqualFlow(; name::Symbol=Symbol("1"), v::Int=0) = new(name, [], [], Ref(v))
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
    Port(node::Junction, index) = new(node, index)
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
    name::Symbol
    nodes::Vector{T} where T <: AbstractNode
    bonds::Vector{Bond}
end
BondGraph(name::Symbol=:BG) = BondGraph(:BG, name, AbstractNode[], Bond[])

# New component
function new(type,name::Symbol=type;library=standard_library)
    d = library[type]
    p = collect(keys(d[:parameters]))
    if haskey(d,:state_vars)
        x = collect(keys(d[:state_vars]))
    else
        x = Vector{Num}([])
    end
    Component(type, name; numports=d[:numports], parameters=p, state_vars=x, equations=d[:equations])
end

# All models
Model = Union{AbstractNode,BondGraph}

# Component type
type(n) = n.type
type(n::EqualEffort) = Symbol("0")
type(n::EqualFlow) = Symbol("1")

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Ports
freeports(n::AbstractNode) = n.freeports
numports(n::AbstractNode) = length(n.freeports)
updateport!(n::AbstractNode, idx::Int) = freeports(n)[idx] = !freeports(n)[idx]

nextfreeport(n::AbstractNode) = findfirst(freeports(n))
function nextfreeport(j::Junction)
    push!(j.freeports,true)
    numports(j)
end
function nextfreeport(j::EqualFlow)
    push!(j.weights,0)
    push!(j.freeports,true)
    numports(j)
end

nextsrcport(n::AbstractNode) = nextfreeport(n)
function nextsrcport(n::EqualFlow)
    i = nextfreeport(n)
    n.weights[i] = -1
    i
end

nextdstport(n::AbstractNode) = nextfreeport(n)
function nextdstport(n::EqualFlow)
    i = nextfreeport(n)
    n.weights[i] = 1
    i
end

# Nodes in Bonds
srcnode(b::Bond) = b.src.node
dstnode(b::Bond) = b.dst.node
in(n::AbstractNode, b::Bond) = n == srcnode(b) || n == dstnode(b)
Base.iterate(b::Bond) = (b.src,true)
function Base.iterate(b::Bond,state)
    if state
        return (b.dst,false)
    else
        return nothing
    end
end

# Parameters
params(n::Component) = n.parameters
params(j::Junction) = Vector{Num}([])

# Components
components(bg::BondGraph) = bg.nodes

# State variables
state_vars(n::Component) = n.state_vars
state_vars(j::Junction) = Vector{Num}([])
function state_vars(m::BondGraph)
    i = 0
    states = Vector{Tuple{Model,Num}}([])
    for c in components(m), x in state_vars(c)
        push!(states,(c,x))
    end
    @variables x[1:length(states)](t)
    return OrderedDict(x[i] => y for (i,y) in enumerate(states))
end

# Set parameter values
function set_param!(n::Component,var,val)
    any(isequal.(params(n),var)) || error("Component does not have parameter.")
    default_value(n)[var] = val
end

# Set initial conditions
function set_initial_value!(n::Component,var,val)
    any(isequal.(state_vars(n),var)) || error("Component does not have state variable.")
    default_value(n)[var] = val
end

# Get default values
default_value(n::Component) = n.default
default_value(n::Component,v::Num) = default_value(n::Component)[v]

# Bonds
bonds(m::BondGraph) = m.bonds

# I/O
show(io::IO, node::AbstractNode) = print(io, "$(type(node)):$(node.name)")
#show(io::IO, node::Junction) = print(io, "$(type(node))")
show(io::IO, port::Port) = print(io, "Port $(port.node) ($(port.index))")
show(io::IO, b::Bond) = print(io, "Bond $(srcnode(b)) ⇀ $(dstnode(b))")
show(io::IO, bg::BondGraph) = print(io, "BondGraph $(type(bg)):$(bg.name) ($(lg.nv(bg)) Nodes, $(lg.ne(bg)) Bonds)")
