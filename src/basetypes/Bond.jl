# PORT
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


# BOND
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

