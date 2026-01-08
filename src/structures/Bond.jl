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
function Bond(srccomp::Component, dstcomp::Component)
    hasfreeport(srccomp) || error("$srccomp has no free ports")
    hasfreeport(dstcomp) || error("$dstcomp has no free ports")
    srcport = nextfreeport(srccomp)
    dstport = nextfreeport(dstcomp)
    Bond(srcport, dstport)
end

ports(b::Bond) = b.src, b.dst
vertices(b::Bond) = parent(b.src), parent(b.dst)

function Base.show(io::IO, b::Bond)
    src, dst = vertices(b)
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
in(n::AbstractNode, b::Bond) = n === srcnode(b) || n === dstnode(b)

iterate(b::Bond) = (b.src, true)
iterate(b::Bond, state) = state ? (b.dst, false) : nothing

# function show(io::IO, b::Bond)
#     print(io, "Bond $(b.src[1])[$(b.src[2])] ⇀ $(b.dst[1])[$(b.dst[2])]")
# end
