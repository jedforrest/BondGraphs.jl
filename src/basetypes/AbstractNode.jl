abstract type AbstractNode end

@enum PortConnection In = -1 Free = 0 Out = 1

# COMPONENT
struct Component{N} <: AbstractNode
    type::Symbol
    name::Symbol
    ports::MVector{N,PortConnection}
    vertex::RefValue{Int}
    parameters::Vector{Num}
    state_vars::Vector{Num}
    equations::Vector{Equation}
    function Component{N}(t, n, v::Int, p::Vector, x::Vector, eq::Vector) where {N}
        portconnections = MVector{N,PortConnection}(repeat([Free], N))
        new(Symbol(t), Symbol(n), portconnections, Ref(v), p, x, eq)
    end
end
# function Component(type, name=type;
#     numports::Int=1, vertex::Int=0, parameters::Vector=[], state_vars::Vector=[],
#     equations::Vector=[])
#     Component{numports}(type, name, numports, vertex, parameters, state_vars, default, equations)
# end

function Component(type, name = type; library = standard_library,
    numports::Int = 1, vertex::Int = 0, parameters::Vector = Num[], state_vars::Vector = Num[],
    equations::Vector = Equation[])

    if !isnothing(library) && type in keys(library)
        d = library[type]
        parameters = collect(keys(d[:parameters]))
        state_vars = if haskey(d, :state_vars)
            collect(keys(d[:state_vars]))
        else
            Num[]
        end
        N = d[:numports]
        Component{N}(type, name, vertex, parameters, state_vars, d[:equations])
    else
        Component{numports}(type, name, vertex, parameters, state_vars, equations)
    end
end


# JUNCTION
abstract type Junction <: AbstractNode end

struct EqualEffort <: Junction
    name::Symbol
    ports::Vector{PortConnection}
    vertex::RefValue{Int}
    EqualEffort(; name = :𝟎, v::Int = 0) = new(Symbol(name), PortConnection[], Ref(v))
end

struct EqualFlow <: Junction
    name::Symbol
    ports::Vector{PortConnection}
    # weights::Vector{Int}
    vertex::RefValue{Int}
    EqualFlow(; name = :𝟏, v::Int = 0) = new(Symbol(name), PortConnection[], Ref(v))
end


# PROPERTIES
# Type
type(n::AbstractNode) = n.type
type(j::Junction) = j.name

# Name
name(n::AbstractNode) = n.name

# Ports
# freeports(n::AbstractNode) = n.freeports
portconnections(n::AbstractNode) = string.(n.ports)
portweights(n::AbstractNode) = Int.(n.ports)
numports(n::AbstractNode) = length(n.ports)
updateport!(n::AbstractNode, idx::Int, pc::PortConnection) = n.ports[idx] = pc

nextfreeport(n::AbstractNode) = findfirst(portconnections(n) .== "Free")
function nextfreeport(j::Junction)
    push!(j.ports, Free)
    numports(j)
end
# function nextfreeport(j::EqualFlow)
#     push!(j.weights,0)
#     push!(j.freeports,true)
#     numports(j)
# end

# nextsrcport(n::AbstractNode) = nextfreeport(n)
# function nextsrcport(n::EqualFlow)
#     i = nextfreeport(n)
#     n.weights[i] = -1
#     i
# end

# nextdstport(n::AbstractNode) = nextfreeport(n)
# function nextdstport(n::EqualFlow)
#     i = nextfreeport(n)
#     n.weights[i] = 1
#     i
# end

# Weights (used when generating equations)
# function bondweights(bg::BondGraph, n::AbstractNode)
#     [nbr in g.inneighbors(bg, n) ? -1 : 1 for nbr in g.all_neighbors(bg, n)]
# end

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Parameters
params(n::AbstractNode) = n.parameters
params(::Junction) = Num[]

# State variables
state_vars(n::AbstractNode) = n.state_vars
state_vars(::Junction) = Num[]

# Equations
equations(n::AbstractNode) = n.equations
equations(::Junction) = Equation[]

# Set parameter values
# function set_param!(n::Component,var,val)
#     any(isequal.(params(n),var)) || error("Component does not have parameter.")
#     default_value(n)[var] = val
# end

# # Set initial conditions
# function set_initial_value!(n::Component,var,val)
#     any(isequal.(state_vars(n),var)) || error("Component does not have state variable.")
#     default_value(n)[var] = val
# end

# # Get default values
# default_value(n::Component) = n.default
# default_value(n::Component,v::Num) = default_value(n::Component)[v]
# default_value(j::Junction) = Dict{Num,Any}()


# BASE FUNCTIONS
# This definition will need to expand when equations etc. are added
==(n1::AbstractNode, n2::AbstractNode) = type(n1) == type(n2) && n1.name == n2.name

show(io::IO, node::AbstractNode) = print(io, "$(node.type):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.name)")
