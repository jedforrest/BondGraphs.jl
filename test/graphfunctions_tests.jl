@testset "Basic Graph Functionality" begin
    c1 = Component(:C)
    c2 = Component(:C, "newC", 1)
    j = Junction(:J)

    b = Bond(c1, j)
    @test src(b) == c1
    @test dst(b) == j

    bg = BondGraph()
    @testset "BondGraph Properties" begin
        @test bg.metamodel == :BG
        @test isempty(bg.nodes)
        @test eltype(bg) == AbstractNode
        @test edgetype(bg) == LightGraphs.AbstractSimpleEdge{Integer}
        @test is_directed(bg)
    end

    @testset "Adding and removing elements" begin
        @test add_vertex!(bg, c1)
        @test add_edge!(bg, b)

        add_vertex!(bg, c2)
        add_vertex!(bg, j)

        @test ne(bg) == 1
        @test has_edge(bg, b)
        @test has_edge(bg, c1, j)
        @test !has_edge(bg, c1, c2)

        @test nv(bg) == 3
        @test has_vertex(bg, j)

        @test inneighbors(bg, c1) == []
        @test outneighbors(bg, c1) == [j]

        @test rem_edge!(bg, b)
        @test ne(bg) == 0
        @test rem_vertex!(bg, c2)
        @test nv(bg) == 2
    end

    @testset "LightGraph Extra Functions" begin

    end
end






