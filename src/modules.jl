function expose(bg::BondGraph, ports::Vector{SourceSensor})
    N = length(ports)
    BondGraphNode(bg, :BG, name(bg), fill(true,N), Ref(0))
end