@testset "Library Functions" begin
    # Standard component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :C)
    # Biochemical component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :Re)
    # Non-existent component
    @test !haskey(BondGraphs.DEFAULT_LIBRARY, :A)

    lib = Dict(:A => Dict())
    addlibrary!(lib)
    # Component now exists
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :A)

    # Delete fake component for later tests
    delete!(BondGraphs.DEFAULT_LIBRARY, :A)
end

@testset "Graph Node Colour Selection" begin
    nodes = [Component(:C), Component(:Re), Component(:SS), EqualFlow()]
    @test BondGraphs.nodecolours(nodes) == [1, 2, 3, :lightgray]
end

@testset "Plotting" begin
    # See conversion_tests.jl
    rn = @reaction_network MM_reversible begin
        (1, 1), E + S <--> C
        (1, 1), C <--> E + P
    end
    bg = BondGraph(rn; chemostats=["S", "P"])
    
    rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), bg)
    attributes = getfield(rec[1], 1)

    @test attributes[:curves] == false
    @test attributes[:title] == :MM_reversible
    @test attributes[:nodeshape] == :rect
end

@testset "Latexify" begin
    bg = RLC()
    eq = equations(bg)
    ltx = repr("text/latex", eq)

    @test ltx == "\\begin{align}\n\\frac{dC_{+}q(t)}{dt} =& C_{+p1_{+}F}\\left( t \\right) \\\\\n\\frac{dI_{+}p(t)}{dt} =& R_{+p1_{+}E}\\left( t \\right) \\\\\n0 =&  - R_{+}R R_{+p1_{+}F}\\left( t \\right) + R_{+p1_{+}E}\\left( t \\right)\n\\end{align}\n"
end