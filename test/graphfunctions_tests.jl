@testset "BondGraph Properties" begin
    bg = BondGraph()
    @test bg.metamodel == :BG
    @test isempty(bg.nodes)
    @test eltype(bg) == AbstractNode
    @test edgetype(bg) == LightGraphs.AbstractSimpleEdge{Integer}
    @test is_directed(bg)
end

@testset "Adding and removing elements" begin
    c1 = Component(:C)
    c2 = Component(:C, "newC", 1)
    j = Junction(:J)

    b = Bond(c1, j)
    @test src(b) == c1
    @test dst(b) == j

    bg = BondGraph()

    @test add_vertex!(bg, c1)
    @test add_edge!(bg, b)

    # adding the same edge twice should fail
    @test !add_edge!(bg, b)

    add_vertex!(bg, c2)
    add_vertex!(bg, j)

    @test ne(bg) == 1
    @test has_edge(bg, b)
    @test has_edge(bg, c1, j)
    @test !has_edge(bg, c1, c2)

    @test nv(bg) == 3
    @test has_vertex(bg, j)

    @test inneighbors(bg, c1) == []
    @test outneighbors(bg, c1) == [3]

    @test rem_edge!(bg, b)
    @test ne(bg) == 0
    @test rem_vertex!(bg, c2)
    @test nv(bg) == 2
end

@testset "LightGraph Extra Functions" begin
    c1 = Component(:C)
    c2 = Component(:R)
    c3 = Component(:I)
    j = Junction(:J)

    b1 = Bond(c1, j)
    b2 = Bond(j, c2)
    b3 = Bond(j, c3)

    bg = BondGraph()

    add_vertex!(bg, c1)
    add_vertex!(bg, c2)
    add_vertex!(bg, c3)
    add_vertex!(bg, j)

    add_edge!(bg, b1)
    add_edge!(bg, b2)
    add_edge!(bg, b3)

    # Testing on a selection of core LG functions
    @test Δ(bg) == 3
    @test density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]
end