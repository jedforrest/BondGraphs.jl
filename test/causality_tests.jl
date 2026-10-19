bg = RCI()

B = incidence_matrix(bg)
@test B == [
    0  1  0  0
    0  0  1  0
    0  0  0  1
    -1 0  0  0
    1 -1 -1 -1
]
@test all(b.causality == :none for b in bonds(bg))

# assign_causality(bg)

# reset causal strokes
for b in bonds(bg)
    b.causality = :none
end

# 1. Assign causalities to Se and Sf components, propagate
sources = BondGraphs.filterbytype(BondGraphs.SourceElement, components(bg))
S = sources[]

getindex(bg, :v)


if S isa EffortSource

elseif S isa FlowSource

end

idx = code_for(bg.graph, name(S))
T = bfs_tree(bg, idx)  # breadth first search

e1 = first(edges(T))
e1.dst
_get_bond(bg, edge) = bg[label_for(bg.graph, edge.src), label_for(bg.graph, edge.dst)]
b = _get_bond(bg, e1)
src = bg[label_for(bg.graph, e1.src)]
dst = bg[label_for(bg.graph, e1.dst)]

function _apply_causality_rule(elem::BondGraphs.AbstractElement)
    if elem isa EffortSource
        return :dst
    elseif elem isa FlowSource
        return :src
    elseif elem isa StaticStorageElement
        # TODO CONTINUE FROM HERE
        # need easy way to determine causality systematically
    else
        @warn "unknown causality rule"
    end
end

_apply_causality_rule(src)


# same graph, but different edge order
for edge in edges(T)
    @show b = _get_bond(bg, edge)

end

# 2. Assign integral causalities to C and I components, propagate

# 3. Assign arbitrary causality to remaining bonds, propagate


@test all(assign_causality(bg) == [:dst, :src, :src, :dst])

C = causality_matrix(bg)
@test C == [
    0  1  0  0
    0  0  1  0
    0  0  0  1
    -1 0  0  0
    1 -1 -1 -1
]
