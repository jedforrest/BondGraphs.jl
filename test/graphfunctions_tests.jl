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
    c2 = Component(:C, name="newC", maxports=1)
    j = Junction(:J)

    b = Bond(c1, j)
    @test src(b) == vertex(c1)
    @test dst(b) == vertex(j)

    bg = BondGraph()

    @test add_vertex!(bg, c1)
    @test add_edge!(bg, b)

    # adding the same edge twice should fail
    @test !add_edge!(bg, b)

    add_vertex!(bg, c2)
    add_vertex!(bg, j)

    @test ne(bg) == 1
    @test has_edge(bg, b)
    @test has_edge(bg, vertex(c1), vertex(j))
    @test !has_edge(bg, vertex(c1), vertex(c2))

    @test nv(bg) == 3
    @test has_vertex(bg, j)

    @test inneighbors(bg, vertex(c1)) == []
    @test outneighbors(bg, vertex(c1)) == [3]

    @test rem_edge!(bg, b)
    @test ne(bg) == 0
    @test rem_vertex!(bg, c2)
    @test nv(bg) == 2
end

@testset "Printing" begin
    C = Component(:C)
    SS = Component(:SS, name="Source")
    J0 = Junction(:J0)
    b1 = Bond(C,J0)
    b2 = Bond(J0,SS)
    bg = BondGraph(name="newbg")

    # repr returns the output of the 'show' function
    @test repr(C) == "C:C"
    @test repr(SS) == "SS:Source"
    @test repr(b1) == "Bond C:C ⇀ J0"
    @test repr(b2) == "Bond J0 ⇀ SS:Source"
    @test repr(bg) == "BondGraph BG:newbg (0 Nodes, 0 Bonds)"

    add_vertex!(bg, C)
    add_vertex!(bg, SS)
    add_vertex!(bg, J0)
    add_edge!(bg, b1)
    add_edge!(bg, b2)
    @test repr(bg) == "BondGraph BG:newbg (3 Nodes, 2 Bonds)"
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

    # Testing on a selection of common LG functions
    @test Δ(bg) == 3
    @test density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]
end