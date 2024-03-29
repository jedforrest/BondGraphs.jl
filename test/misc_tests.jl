@testset "Library Functions" begin
    # Standard component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :C)
    # Biochemical component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :Re)
    # Non-existent component
    @test !haskey(BondGraphs.DEFAULT_LIBRARY, :A)

    A = Dict(
        :description => "Some description",
        :numports => 5,
    )

    lib = Dict(:A => A)
    addlibrary!(lib)
    # Component now exists
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :A)
    @test numports(Component(:A)) == 5

    # Delete fake components for later tests
    delete!(BondGraphs.DEFAULT_LIBRARY, :A)
end

@testset "Graph Attributes" begin
    nodes = [
        Component(:C),
        Component(:Re, "R1"),
        SourceSensor(),
        EqualFlow(),
        EqualEffort(),
        BondGraphNode(BondGraph("testBG"))
    ]

    @test BondGraphs.nodecolour.(nodes) == [:lightblue, :coral, :lightgreen, :white, :white, :white]
    @test BondGraphs.nodelabel.(nodes) == ["C", "R1", "SS", "1", "0", "testBG"]
end

@testset "Plotting" begin
    # See catalyst_tests.jl
    rn = @reaction_network MM_reversible begin
        (1, 1), E + S <--> C
        (1, 1), C <--> E + P
    end
    bg = BondGraph(rn; chemostats=["S", "P"])

    rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bg)
    attributes = getfield(rec[1], 1)

    @test attributes[:curves] == false
    @test attributes[:title] == "MM_reversible"
    @test attributes[:nodeshape] == :circle
end

# @testset "Latexify" begin
#     bg = RLC()
#     eq = equations(bg)
#     ltx = repr("text/latex", eq)

#     @test ltx == "\\begin{align}\n\\frac{dC_{+}q(t)}{dt} =& C_{+p1_{+}F}\\left( t \\right) \\\\\n\\frac{dI_{+}p(t)}{dt} =& R_{+p1_{+}E}\\left( t \\right)\n\\end{align}\n"
# end
