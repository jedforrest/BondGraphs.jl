function expose(bg::BondGraph, ports::Vector{SourceSensor})
    N = length(ports)
    BondGraphNode(bg, :BG, name(bg), ports, fill(true,N), Ref(0))
end

exposed(bgn) = bgn.exposed