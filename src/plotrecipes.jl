# @recipe function f(bg::BondGraph)
#     # curves --> false

#     # GraphRecipes.graphplot(bg, curves = false, title = "test")
#     # RecipesBase.recipetype(:graphplot, bg)
#     # graphplot(bg)
#     # plt
#     ()
# end

@recipe function f(bg::BondGraph)
    markershape --> :auto        # if markershape is unset, make it :auto
    xrotation --> 45           # if xrotation is unset, make it 45
    zrotation --> 90           # if zrotation is unset, make it 90
    # rand(10, 10)
    bg
end