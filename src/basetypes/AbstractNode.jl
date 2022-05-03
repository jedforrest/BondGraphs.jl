abstract type AbstractNode end

# COMPONENT
struct Component{N} <: AbstractNode
    type::Symbol
    name::Symbol
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    variables::Dict{Symbol,Dict{Num,Any}}
    equations::Vector{Equation}
    function Component{N}(t, n, vx, vars, eq) where {N}
        new(Symbol(t), Symbol(n), ones(MVector{N,Bool}), Ref(vx), vars, eq)
    end
end

function Component(type, name=type; vertex::Int=0, library=BondGraphs.DEFAULT_LIBRARY,
    comp_dict=_get_comp_default(library, type),
    numports::Int=_get_comp_default(comp_dict, :numports, 1),
    vars=_get_comp_default(comp_dict, :variables),
    equations=_get_comp_default(comp_dict, :equations, Equation[]),
    kwargs...)

    # add default empty dicts to variables dict
    vars_empty = Dict(:parameters=>Dict(), :globals=>Dict(), :states=>Dict(), :controls=>Dict())
    vars = deepcopy(merge(vars_empty, vars))

    # Actual construction of the component
    comp = Component{numports}(type, name, vertex, vars, equations)

    # kwargs are used to set default variable values
    for (k,v) in kwargs
        setproperty!(comp, k, v)
    end

    comp
end

_get_comp_default(D, key, default=Dict()) = haskey(D, key) ? D[key] : default

# Source-sensor
struct SourceSensor <: AbstractNode
    name::Symbol
    freeports::MVector{1,Bool}
    vertex::RefValue{Int}
    SourceSensor(; name=:SS, v::Int=0) = new(Symbol(name), ones(MVector{1,Bool}), Ref(v))
end


# JUNCTION
abstract type Junction <: AbstractNode end

struct EqualEffort <: Junction
    name::Symbol
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    EqualEffort(; name=Symbol("0"), v::Int=0) = new(Symbol(name), [true], [0], Ref(v))
end

struct EqualFlow <: Junction
    name::Symbol
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    EqualFlow(; name=Symbol("1"), v::Int=0) = new(Symbol(name), [true], [0], Ref(v))
end


# PROPERTIES
# Type
type(n::AbstractNode) = n.type
type(j::Junction) = typeof(j)
type(::SourceSensor) = :SS

# Name
name(n::AbstractNode) = n.name

# Ports
freeports(n::AbstractNode) = n.freeports
numports(n::AbstractNode) = length(n.freeports)
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
parameters(::AbstractNode) = ()
parameters(n::Component) = n.variables[:parameters]

# Globals
globals(::AbstractNode) = ()
globals(n::Component) = n.variables[:globals]

# State variables
states(::AbstractNode) = ()
states(n::Component) = n.variables[:states]

# Control variables
controls(::AbstractNode) = ()
controls(n::Component) = n.variables[:controls]

# Equations
equations(::AbstractNode) = Equation[]
equations(n::Component) = n.equations

# Defaults
# defaults(n::AbstractNode) = merge(n.variables.parameters, n.variables.globals,
#     n.variables.states, n.variables.controls)
# defaults(::Junction) = Dict{Num,Any}()
# defaults(::SourceSensor) = Dict{Num,Any}()

# Variables
all_variables(::AbstractNode) = ()
function all_variables(n::Component)
    # use getfield here so that this can be used by getproperty
    merge(values(getfield(n, :variables))...)
end

# function get_default(n::AbstractNode, var)
#     default_dict, _var = _find_var(n, var)
#     isnothing(_var) && error("Could not find '$var' in node $n")
#     default_dict[_var]
# end
# function set_default!(n::AbstractNode, var, val)
#     default_dict, _var = _find_var(n, var)
#     isnothing(_var) && error("Could not find '$var' in node $n")
#     if default_dict isa Dict{Num,Function} && val isa Number
#         # control var constant replace by constant function
#         default_dict[_var] = (t -> val)
#     else
#         default_dict[_var] = val
#     end
# end

# function _find_var(n::AbstractNode, var)
#     # Try both parameters and state/control vars
#     p, = @parameters $var
#     _, x = @variables t, $var(t)
#     if string(p) in string.(parameters(n))
#         return n.parameters, p
#     elseif string(p) in string.(globals(n))
#         return n.globals, GlobalScope(p)
#     elseif string(x) in string.(states(n))
#         return n.states, x
#     elseif string(x) in string.(controls(n))
#         return n.controls, x
#     else
#         return nothing, nothing
#     end
# end

# BASE FUNCTIONS
# This definition will need to expand when equations etc. are added
==(n1::AbstractNode, n2::AbstractNode) = type(n1) == type(n2) && n1.name == n2.name

show(io::IO, node::AbstractNode) = print(io, "$(type(node)):$(node.name)")
show(io::IO, node::Junction) = print(io, "$(node.name)")

# Easier referencing systems using a.b notation
function getproperty(n::Component, sym::Symbol)
    p, = @parameters $sym
    _, x = @variables t, $sym(t)
    all_vars = all_variables(n)
    if p in keys(all_vars)
        return all_vars[p]
    elseif x in keys(all_vars)
        return all_vars[x]
    else
        getfield(n, sym)
    end
end

function setproperty!(n::Component, sym::Symbol, val)
    p, = @parameters $sym
    _, x = @variables t, $sym(t)

    for (_, vars) in getfield(n, :variables)
        if p in keys(vars)
            return vars[p] = val
        elseif x in keys(vars)
            return vars[x] = val
        end
    end
    setfield!(n, sym, val)
end
