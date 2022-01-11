abstract type AbstractNode end

# COMPONENT
struct Component{N} <: AbstractNode
    type::Symbol
    name::Symbol
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    parameters::Vector{Num} # store as dict? Num => Real
    states::Vector{Num}
    equations::Vector{Equation}
    function Component{N}(t, n, v::Int, p::Vector, x::Vector, eq::Vector) where {N}
        new(Symbol(t), Symbol(n), ones(MVector{N,Bool}), Ref(v), p, x, eq)
    end
end

function Component(type, name = type; library = BondGraphs.DEFAULT_LIBRARY[],
    numports::Int = 1, vertex::Int = 0, parameters::Vector = Num[], states::Vector = Num[],
    equations::Vector = Equation[])

    if !isnothing(library) && type in keys(library)
        d = library[type]
        parameters = collect(keys(d[:parameters]))
        states = if haskey(d, :states)
            collect(keys(d[:states]))
        else
            Num[]
        end
        numports = d[:numports]
        equations = d[:equations]
    end
    Component{numports}(type, name, vertex, parameters, states, equations)
end


# JUNCTION
abstract type Junction <: AbstractNode end

struct EqualEffort <: Junction
    name::Symbol
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    EqualEffort(; name = Symbol("0"), v::Int = 0) = new(Symbol(name), [true], [0], Ref(v))
end

struct EqualFlow <: Junction
    name::Symbol
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    EqualFlow(; name = Symbol("1"), v::Int = 0) = new(Symbol(name), [true], [0], Ref(v))
end


# PROPERTIES
# Type
type(n::AbstractNode) = n.type
type(j::Junction) = typeof(j)

# Name
name(n::AbstractNode) = n.name

# Ports
freeports(n::AbstractNode) = n.freeports
numports(n::AbstractNode) = length(n.freeports)
# numports(::Junction) = Inf
updateport!(n::AbstractNode, idx::Int) = freeports(n)[idx] = !freeports(n)[idx]

# Weights
weights(j::Junction) = j.weights
set_weight!(j::Junction, idx::Int, w::Int) = j.weights[idx] = w

nextfreeport(n::AbstractNode) = findfirst(freeports(n))
function nextfreeport(j::Junction)
    freeport = findfirst(freeports(j))
    if isnothing(freeport)
        push!(j.freeports, true)
        push!(j.weights, 0)
        return length(j.freeports)
    else
        return freeport
    end
end

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Parameters
parameters(n::AbstractNode) = n.parameters
parameters(::Junction) = Num[]

# State variables
states(n::AbstractNode) = n.states
states(::Junction) = Num[]

# Equations
equations(n::AbstractNode) = n.equations
equations(::Junction) = Equation[]

# Set parameter values
# function set_param!(n::Component,var,val)
#     any(isequal.(parameters(n),var)) || error("Component does not have parameter.")
#     default_value(n)[var] = val
# end

# # Set initial conditions
# function set_initial_value!(n::Component,var,val)
#     any(isequal.(states(n),var)) || error("Component does not have state variable.")
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
