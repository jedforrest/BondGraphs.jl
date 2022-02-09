@recipe function f(bg::BondGraph)
    # Default attributes
    nodecolor --> nodecolours(bg.nodes)
    title --> bg.name
    names --> bg.nodes

    # Forced attributes
    curves := false
    nodeshape := :rect

    g.SimpleDiGraph(g.adjacency_matrix(bg))
end

# @recipe function f(bg::BondGraph)
#     curves := false
#     RecipesBase.recipetype(:graphplot, bg)
# end

function nodecolours(nodes)
    colours = []
    for n in nodes
        type_n = type(n)
        # colour i is the i-th default plotting colour
        c = if type_n in [:Se, :Sf, :SS]
            3
        elseif n isa Component{1}
            1
        elseif n isa Component # more than 1 port
            2
        elseif n isa Junction
            :lightgray
        end
        push!(colours, c)
    end
    colours
end
