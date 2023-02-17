abstract type AbstractNode end

"""
    Component{N} <: AbstractNode

`N`-port bond graph component.
"""
struct Component{N} <: AbstractNode
    type::AbstractString
    name::AbstractString
    freeports::MVector{N,Bool}
    vertex::RefValue{Int}
    variables::Dict{Symbol,Dict{Num,Any}}
    equations::Vector{Equation}
    function Component{N}(t, n, vx, vars, eq) where {N}
        new(string(t), string(n), ones(MVector{N,Bool}), Ref(vx), vars, eq)
    end
end

"""
    Component(type, name=type)
    Component(type, name=type; library=BondGraphs.DEFAULT_LIBRARY, <keyword arguments>)

Construct a Component of a defined (bondgraph) type ∈ {R, C, I, Se, Sf, TF, Ce, Re, SCe}.
Properties and equations of available componets are defined in `BondGraphs.DEFAULT_LIBRARY`.
"""
function Component(type, name=type;
    vertex::Int=0,
    library=BondGraphs.DEFAULT_LIBRARY,
    comp_dict=_get_comp_default(library, type),
    numports::Int=_get_comp_default(comp_dict, :numports, 1),
    vars=_get_comp_default(comp_dict, :variables),
    equations=_get_comp_default(comp_dict, :equations, Equation[]),
    kwargs...)

    # add default empty dicts to variables dict
    vars_empty = Dict(:parameters => Dict(), :globals => Dict(), :states => Dict(), :controls => Dict())
    vars = deepcopy(merge(vars_empty, vars))

    # Actual construction of the component
    comp = Component{numports}(type, name, vertex, vars, equations)

    # kwargs are used to set default variable values
    for (k, v) in kwargs
        setproperty!(comp, k, v)
    end

    comp
end

_get_comp_default(D, key, default=Dict()) = haskey(D, key) ? D[key] : default

"""
    SourceSensor <: AbstractNode

Special component type that acts as a source of both effort and flow
"""
struct SourceSensor <: AbstractNode
    name::AbstractString
    freeports::MVector{1,Bool}
    vertex::RefValue{Int}
    function SourceSensor(; name="SS", v::Int=0)
        new(string(name), ones(MVector{1,Bool}), Ref(v))
    end
end


# JUNCTION
abstract type Junction <: AbstractNode end

"""
    EqualEffort <: Junction

Efforts are all equal, flows sum to zero (aka 0-junction).
"""
struct EqualEffort <: Junction
    name::AbstractString
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    function EqualEffort(; name="𝟎", v::Int=0)
        new(string(name), [true], [0], Ref(v))
    end
end

"""
    EqualFlow <: Junction

Flows are all equal, efforts sum to zero (aka 1-junction).
"""
struct EqualFlow <: Junction
    name::AbstractString
    freeports::Vector{Bool}
    weights::Vector{Int}
    vertex::RefValue{Int}
    function EqualFlow(; name="𝟏", v::Int=0)
        new(string(name), [true], [0], Ref(v))
    end
end


# PROPERTIES
# Type
type(n::AbstractNode) = n.type
type(j::Junction) = typeof(j)
type(::SourceSensor) = "SS"

# Name
name(n::AbstractNode) = n.name
name(n::Junction) = vertex(n) == 0 ? n.name : "$(n.name)_$(vertex(n))"

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
parameters(::AbstractNode) = Dict()
parameters(n::Component) = n.variables[:parameters]

# Globals
globals(::AbstractNode) = Dict()
globals(n::Component) = n.variables[:globals]

# State variables
states(::AbstractNode) = Dict()
states(n::Component) = n.variables[:states]

# Control variables
controls(::AbstractNode) = Dict()
controls(n::Component) = n.variables[:controls]

# Equations
equations(::AbstractNode) = Equation[]
equations(n::Component) = n.equations

# Variables
all_variables(::AbstractNode) = ()
function all_variables(n::Component)
    # use getfield here so that this can be used by getproperty
    merge(values(getfield(n, :variables))...)
end

# BASE FUNCTIONS
# This definition will need to expand when equations etc. are added
==(n1::AbstractNode, n2::AbstractNode) = type(n1) == type(n2) && name(n1) == name(n2)

show(io::IO, node::AbstractNode) = print(io, "$(type(node)):$(name(node))")
show(io::IO, node::Junction) = print(io, name(node))

# Easier referencing systems using a.b notation
# TODO: rearrange so that getfield() is checked first with isdefined()
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

# TODO: this can overrite global variables unintentionally by creating a new local copy
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
