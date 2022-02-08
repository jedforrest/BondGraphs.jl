@testset "Library Functions" begin
    # Standard component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :C)
    # Biochemical component
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :Re)
    # Non-existing component
    @test !haskey(BondGraphs.DEFAULT_LIBRARY, :A)

    lib = Dict(:A => Dict())
    addlibrary!(lib)
    # Component now exists
    @test haskey(BondGraphs.DEFAULT_LIBRARY, :A)
end

@testset "Graph Node Colour Selection" begin
    nodes = [Component(:C), Component(:Re), Component(:SS), EqualFlow()]
    @test BondGraphs.nodecolours(nodes) == [1, 2, 3, :lightgray]
end