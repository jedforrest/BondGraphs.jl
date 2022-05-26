@recipe function f(bg::BondGraph)
    # Default attributes
    nodecolor --> nodecolours(bg.nodes)
    title --> bg.name
    names --> nodenames(bg.nodes)

    # Forced attributes
    curves := false
    nodeshape := :rect

    GraphRecipes.GraphPlot([bg])
end

function nodecolours(nodes)
    # colour integer i is the i-th default plotting colour
    [
        (if type(n) in ["Se", "Sf", "SS", "SCe"] # Source node types
            3
        elseif n isa Component{1}
            1
        elseif n isa Component # more than 1 port
            2
        elseif n isa Junction
            :lightgray
        else
            :white # blank
        end)
        for n in nodes
    ]
end

function nodenames(nodes)
    [
        (if n isa EqualEffort
            "0"
        elseif n isa EqualFlow
            "1"
        else
            "$(type(n)):$(name(n))"
        end)
        for n in nodes
    ]
end
