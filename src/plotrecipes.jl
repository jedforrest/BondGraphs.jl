# Plots.jl backend
@recipe function plot(bg::BondGraph, showtype=false)
    allnodes = nodes(bg)

    # Default attributes
    nodecolor --> nodecolour.(allnodes)
    title --> name(bg)
    names --> nodelabel.(allnodes, showtype)

    # Forced attributes
    curves := false
    nodeshape := :circle
    strokewidth := 0

    GraphRecipes.GraphPlot([bg])
end

@recipe function plot(bgn::BondGraphNode, showtype=false)
    bg = bgn.bondgraph
    allnodes = nodes(bg)

    # Default attributes
    nodecolor --> nodecolour.(allnodes)
    title --> name(bg)
    names --> nodelabel.(allnodes, showtype)

    # Forced attributes
    curves := false
    nodeshape := :circle
    strokewidth := 0

    GraphRecipes.GraphPlot([bg])
end

function nodecolour(node)
    if type(node) in ["Se", "Sf", "SS", "SCe"] # Source node types
        :lightgreen
    elseif node isa Component{1}
        :lightblue
    elseif node isa Component # more than 1 port
        :coral
    elseif node isa Junction
        :white
    else
        :white
    end
end

function nodelabel(node, showtype=false)
    if node isa EqualEffort
        "0"
    elseif node isa EqualFlow
        "1"
    else
        showtype ? repr(node) : name(node)
    end
end

# Makie/GraphMakie backend
function GraphMakie.graphplot(bg::BondGraph; interactions=true, ongrid=false, showtype=false, kwargs...)
    allnodes = nodes(bg)
    nodesize = 50 .* [n isa Junction ? 0.75 : 1. for n in allnodes]

    fig, ax, p = GraphMakie.graphplot(g.SimpleDiGraph(bg);
        layout=Stress(),
        arrow_shift = :end,
        node_size = nodesize,
        node_color = nodecolour.(allnodes),
        nlabels = nodelabel.(allnodes, showtype),
        nlabels_align=(:center, :center),
        kwargs...
    )

    if ongrid
        rounded_pos = []
        # try rounding to nearest 0, 0.5, or 0.25 successively
        for digits in 0:2
            rounded_pos = [round.(pts; base=2, digits) for pts in p.node_pos[]]
            allunique(rounded_pos) && break
        end
        # layout must be a function with one argument
        newlayout(_) = rounded_pos
        p.layout = newlayout
        GraphMakie.autolimits!(ax)
    end

    # interactions
    GraphMakie.deregister_interaction!(ax, :rectanglezoom)
    if interactions
        GraphMakie.register_interaction!(ax, :nhover, NodeHoverHighlight(p, 1.2))
        GraphMakie.register_interaction!(ax, :ndrag, NodeDrag(p))
    end

    # remove plot decorations and adjust aspect ratio
    GraphMakie.hidedecorations!(ax)
    GraphMakie.hidespines!(ax)
    ax.aspect = GraphMakie.DataAspect()

    fig
end
