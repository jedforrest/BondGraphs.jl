abstract type AbstractNode end

# COMPONENT
struct Component{N} <: AbstractNode
    type::Symbol
    name::Symbol
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    parameters::Vector{Num}
    state_vars::Vector{Num}
    equations::Vector{Equation}
    function Component{N}(t, n, np::Int, v::Int,
         p::Vector, x::Vector, eq::Vector) where N
        new(Symbol(t), Symbol(n), ones(MVector{np,Bool}), Ref(v), p, x, d, eq)
    end
end
# function Component(type, name=type;
#     numports::Int=1, vertex::Int=0, parameters::Vector=[], state_vars::Vector=[],
#     equations::Vector=[])
#     Component{numports}(type, name, numports, vertex, parameters, state_vars, default, equations)
# end

function Component(type, name=type; library=standard_library, 
        numports::Int=1, vertex::Int=0, parameters::Vector=Num[], state_vars::Vector=Num[],
        equations::Vector=Equation[])
    
    if isnothing(library)
        Component{numports}(type, name, numports, vertex, parameters, state_vars, equations)
    else
        d = library[type]
        p = collect(keys(d[:parameters]))
        if haskey(d, :state_vars)
            x = collect(keys(d[:state_vars]))
        else
            x = Num[]
        end
        Component(type, name; numports=d[:numports], parameters=p, state_vars=x, equations=d[:equations])
    end
end


# JUNCTION
abstract type Junction <: AbstractNode end

struct EqualEffort <: Junction
    name::Symbol
    freeports::Vector{Bool}
    vertex::RefValue{Int}
    EqualEffort(; name=:𝟎, v::Int=0) = new(Symbol(name), Bool[], Ref(v))
end

struct EqualFlow <: Junction
    name::Symbol
    freeports::Vector{Bool}
    # weights::Vector{Int}
    vertex::RefValue{Int}
    EqualFlow(; name=:𝟏, v::Int=0) = new(Symbol(name), Bool[], Ref(v))
end


# PROPERTIES
# Type
type(n::AbstractNode) = n.type
type(j::Junction) = j.name

# Name
name(n::AbstractNode) = n.name

# Ports
freeports(n::AbstractNode) = n.freeports
numports(n::AbstractNode) = length(n.freeports)
updateport!(n::AbstractNode, idx::Int) = freeports(n)[idx] = !freeports(n)[idx]

nextfreeport(n::AbstractNode) = findfirst(freeports(n))
function nextfreeport(j::Junction)
    push!(j.freeports, true)
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

# Vertex
vertex(n::AbstractNode) = n.vertex[]
set_vertex!(n::AbstractNode, v::Int) = n.vertex[] = v

# Parameters
params(n::AbstractNode) = n.parameters
params(::Junction) = Num[]

# State variables
state_vars(n::Component) = n.state_vars
state_vars(::Junction) = Num[]

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
