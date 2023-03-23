"""
    Bond(source::AbstractNode, destination::AbstractNode)
    Bond(source::Port, destination::Port)

Connect two bond graph components (or two ports of two components) with a bond. The bond
direction is from `source` to `destination`. If the ports are not specified, the bond will
be created between the next available ports in each component.

In most cases it is better to use [`connect!`](@ref) instead.
"""
struct Bond <: g.AbstractSimpleEdge{Int}
    src::Tuple{AbstractNode, Any}
    dst::Tuple{AbstractNode, Any}
end
function Bond(srcnode::AbstractNode, dstnode::AbstractNode)
    Bond((srcnode, nextfreeport(srcnode)), (dstnode, nextfreeport(dstnode)))
end
# struct Bond <: g.AbstractSimpleEdge{Int}
#     srcport::Port
#     dstport::Port
# end
# function Bond(srcnode::AbstractNode, dstnode::AbstractNode)
#     Bond(Port(srcnode), Port(dstnode))
# end

# Source and Destination
srcnode(b::Bond) = b.src[1]
dstnode(b::Bond) = b.dst[1]

srclabel(b::Bond) = b.src[2]
dstlabel(b::Bond) = b.dst[2]

# Base functions
in(n::AbstractNode, b::Bond) = n === srcnode(b) || n === dstnode(b)

iterate(b::Bond) = (b.src, true)
iterate(b::Bond, state) = state ? (b.dst, false) : nothing

show(io::IO, b::Bond) = print(io, "Bond $(b.src[1])[$(b.src[2])] ⇀ $(b.dst[1])[$(b.dst[2])]")
