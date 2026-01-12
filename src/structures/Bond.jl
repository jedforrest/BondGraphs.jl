"""
    Bond(source::AbstractNode, destination::AbstractNode)
    Bond((source_node, port_label), (destination_node, port_label))

Connect two bond graph components (or two ports of two components) with a bond. The bond
direction is from `source` to `destination`. If the ports are not specified, the bond will
be created between the next available ports in each component.

In most cases it is better to use [`connect!`](@ref) instead.
"""
struct Bond
    src::Port
    dst::Port
    function Bond(srcport::Port, dstport::Port)
        is_connected(srcport) && error("$srcport already connected")
        is_connected(dstport) && error("$dstport already connected")
        connect!(srcport)
        connect!(dstport)
        new(srcport, dstport)
    end
end
function Bond(srccomp::AbstractElement, dstcomp::AbstractElement)
    hasfreeport(srccomp) || error("$srccomp has no free ports")
    hasfreeport(dstcomp) || error("$dstcomp has no free ports")
    srcport = nextfreeport(srccomp)
    dstport = nextfreeport(dstcomp)
    Bond(srcport, dstport)
end

ports(b::Bond) = b.src, b.dst
componentnames(b::Bond) = parentname(b.src), parentname(b.dst)

function Base.show(io::IO, b::Bond)
    src, dst = componentnames(b)
    print(io, "$src ⇀ $dst")
end

# MTK system connector
function connection_equation(b::Bond)
    srcport_sys = system(b.src)
    dstport_sys = system(b.dst)
    ModelingToolkit.connect(srcport_sys, dstport_sys)
end

# # Source and Destination
# srcnode(b::Bond) = b.src[1]
# dstnode(b::Bond) = b.dst[1]

# srclabel(b::Bond) = b.src[2]
# dstlabel(b::Bond) = b.dst[2]

# Base functions
Base.in(n::AbstractElement, b::Bond) = name(n) in componentnames(b)

iterate(b::Bond) = (b.src, true)
iterate(b::Bond, state) = state ? (b.dst, false) : nothing

# src, dst (from Graphs)
# src(b::Bond) = vertex(srcnode(b))
# dst(b::Bond) = vertex(dstnode(b))

# get unique components from a vector of bonds
# FIXME return full components instead of just names
function unique_components(bonds::Vector{Bond})
    return Set(componentnames(b) for b in bonds)
end

disconnect!(b::Bond) = disconnect!.(ports(b))
