abstract type AbstractNode end

struct Component <: AbstractNode
    metamodel::Symbol
    name::String
    maxports::Int
end
Component(metamodel::Symbol, name::String=string(metamodel)) = Component(metamodel, name, 1)

struct Junction <: AbstractNode
    metamodel::Symbol
end

struct Bond <: lg.AbstractSimpleEdge{Integer}
    src::AbstractNode
    dst::AbstractNode
end

struct BondGraph <: lg.AbstractGraph{Integer}
    metamodel::Symbol
    name::String
    nodes::Vector{T} where T <:AbstractNode
    bonds::Vector{Bond}
end
BondGraph() = BondGraph(:BG, "bg", AbstractNode[], Bond[])