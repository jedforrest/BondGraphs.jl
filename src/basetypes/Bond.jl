"""
    Port(node::AbstractNode)
    Port(node::AbstractNode, index::Int)

Create a new Port for `node`. Ports have an index corresponding to the component's variables.

Ports are the `node` elements that are connected by bonds. The port does not technically
exist until this is called, even though a component has a fixed number of assigned ports
when created.

WARNING: connecting a bond to the wrong port may assign values to the wrong variables!
"""
struct Port
    node::AbstractNode
    index::Int
    function Port(node::AbstractNode, index)
        ports = freeports(node)
        any(ports) || error("Node $node has no free ports")
        ports[index] || error("Port $index in node $node is already connected")
        new(node, index)
    end
    Port(node::Junction, index) = new(node, index)
end
Port(node::AbstractNode) = Port(node, nextfreeport(node))


"""
    Bond(source::AbstractNode, destination::AbstractNode)
    Bond(source::Port, destination::Port)

Connect two bond graph components (or two ports of two components) with a bond. The bond
direction is from `source` to `destination`. If the ports are not specified, the bond will
be created between the next available ports in each component.

In most cases it is better to use [`connect!`](@ref) instead.
"""
struct Bond <: g.AbstractSimpleEdge{Int}
    srcport::Port
    dstport::Port
end
function Bond(srcnode::AbstractNode, dstnode::AbstractNode)
    Bond(Port(srcnode), Port(dstnode))
end


# Source and Destination
srcnode(b::Bond) = b.srcport.node
dstnode(b::Bond) = b.dstport.node

# Base functions
in(n::AbstractNode, b::Bond) = n === srcnode(b) || n === dstnode(b)

iterate(b::Bond) = (b.srcport, true)
iterate(b::Bond, state) = state ? (b.dstport, false) : nothing

show(io::IO, port::Port) = print(io, "Port $(port.node) ($(port.index))")
show(io::IO, b::Bond) = print(io, "Bond $(srcnode(b)) ⇀ $(dstnode(b))")
