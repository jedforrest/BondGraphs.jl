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
        srcport.connected = true
        dstport.connected = true
        srcport.weight = OUTBOUND
        dstport.weight = INBOUND
        new(srcport, dstport)
    end
end
# either port or component are accepted as inputs
function Bond(src::Union{Port,AbstractElement}, dst::Union{Port,AbstractElement})
    srcport = src isa Port ? src : nextfreeport(src)
    dstport = dst isa Port ? dst : nextfreeport(dst)
    isnothing(srcport) && error("Component $src has no free ports")
    isnothing(dstport) && error("Component $dst has no free ports")
    Bond(srcport, dstport)
end

ports(b::Bond) = b.src, b.dst
componentnames(b::Bond) = parentname(b.src), parentname(b.dst)

function Base.show(io::IO, b::Bond)
    src, dst = componentnames(b)
    print(io, "$src ⇀ $dst")
end

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

function disconnect!(b::Bond)
    for p in ports(b)
        p.connected = false
        p.weight = DISCONNECTED
    end
end

function connection_equation(b::Bond)
    srcport = system(b.src)
    dstport = system(b.dst)
    ModelingToolkit.connect(srcport, dstport)
end
