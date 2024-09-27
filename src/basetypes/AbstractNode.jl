abstract type AbstractNode end

"""
    Component{N} <: AbstractNode
    Component(type, name=type)
    Component(type, name=type; library=BondGraphs.DEFAULT_LIBRARY, <keyword arguments>)

Construct a Component of a defined (bondgraph) type ∈ {R, C, I, Se, Sf, TF, Ce, Re, SCe}.

Components have a `N` fixed ports when generated. This is usually determined by the bond
graph type. Other properties and equations of available components are defined in
`BondGraphs.DEFAULT_LIBRARY` (see  [`description`](@ref)).
"""
struct Component{N} <: AbstractNode
    type::AbstractString
    name::AbstractString
    ports::OrderedDict{Any,Bool}
    vertex::RefValue{Int}
    variables::Dict{Symbol,Dict{Num,Any}}
    equations::Vector{Equation}
    function Component{N}(t, n, vx, vars, eq) where {N}
        ports = Dict(i => false for i in 1:N)
        new(string(t), string(n), ports, Ref(vx), vars, eq)
    end
end

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

Special component type that acts as a source of both effort and flow. SourceSensors are used
as external ports for [`BondGraphNode`](@ref)s.
"""
struct SourceSensor <: AbstractNode
    name::AbstractString
    ports::Dict{Any,Bool} # might be redundant
    vertex::RefValue{Int}
    function SourceSensor(; name="SS", v::Int=0)
        new(string(name), Dict(1 => false), Ref(v))
    end
end


# JUNCTION
abstract type Junction <: AbstractNode end

"""
    EqualEffort <: Junction

Efforts are all equal, flows sum to zero (0-junction). Has an unlimited number of ports.
"""
struct EqualEffort <: Junction
    name::AbstractString
    ports::Vector{Int} # port number => weight (+1 or -1)
    vertex::RefValue{Int}
    function EqualEffort(; name="𝟎", v::Int=0)
        new(string(name), [0], Ref(v))
    end
end

"""
    EqualFlow <: Junction

Flows are all equal, efforts sum to zero (1-junction). Has an unlimited number of ports.
"""
struct EqualFlow <: Junction
    name::AbstractString
    ports::Vector{Int} # port number => weight (+1 or -1)
    vertex::RefValue{Int}
    function EqualFlow(; name="𝟏", v::Int=0)
        new(string(name), [0], Ref(v))
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
ports(n::AbstractNode) = n.ports
numports(n::AbstractNode) = length(ports(n))
portlabels(n::AbstractNode) = collect(keys(ports(n)))

isconnected(n::AbstractNode, label) = ports(n)[label] != 0 # '!=0' needed for junctions
updateport!(n::AbstractNode, label) = ports(n)[label] = !ports(n)[label]
updateport!(::Junction, ::Int) = nothing # override

port_info(n::AbstractNode) = (n, nextfreeport(n))
port_info(t::Tuple{AbstractNode,Any}) = t
function port_info(t::Tuple{BondGraphNode,String})
    pts = [n for n in nodes(t[1].bondgraph) if n isa SourceSensor]
    for (i,c) in enumerate(pts)
        if (c isa SourceSensor) && (t[2] == c.name)
            return (t[1],i)
        end
    end
    return error("Port $(t[2]) not found.")
end
port_info(t::Tuple{BondGraphNode,Symbol}) = port_info((t[1],string(t[2])))

# ports renamed as ports to make purpose clearer
@deprecate freeports(n::AbstractNode) ports(n::AbstractNode)

# Weights
weights(j::Junction) = j.ports # deprecated
set_weight!(j::Junction, idx::Int, w::Int) = ports(j)[idx] = w

nextfreeport(n::AbstractNode) = findfirst(!, ports(n)) # first 'not' connected port
function nextfreeport(j::Junction)
    index = findfirst(==(0), ports(j))
    if isnothing(index)
        push!(j.ports, 0) # add new empty port
        return numports(j)
    else
        return index
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
