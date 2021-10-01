@testset "BondGraph Properties" begin
    bg = BondGraph()
    @test bg.name == "BG"
    @test isempty(bg.nodes)
    @test eltype(bg) == AbstractNode
    @test edgetype(bg) == LightGraphs.AbstractSimpleEdge{Integer}
    @test is_directed(bg)
end

@testset "Adding and removing elements" begin
    c1 = Component(:C)
    c2 = Component(:C, "newC", numports=1)
    j = Junction(:J)

    b = Bond(c1, j)
    @test src(b) == vertex(c1)
    @test dst(b) == vertex(j)

    bg = BondGraph()

    @test add_vertex!(bg, c1)
    @test add_edge!(bg, c1, j) == b

    add_vertex!(bg, c2)
    add_vertex!(bg, j)

    @test ne(bg) == 1
    @test has_edge(bg, vertex(c1), vertex(j))
    @test !has_edge(bg, vertex(c1), vertex(c2))

    @test nv(bg) == 3
    @test has_vertex(bg, j)

    @test inneighbors(bg, vertex(c1)) == []
    @test outneighbors(bg, vertex(c1)) == [3]

    @test rem_edge!(bg, c1, j) == b
    @test ne(bg) == 0
    @test rem_vertex!(bg, c2)
    @test nv(bg) == 2
end

@testset "BondGraphNode" begin
    bg = BondGraph("RCI")
    bgn = BondGraphNode(bg)

    @test bgn.type == :BG
    @test bgn.name == "RCI"
    @test bgn.freeports == Bool[]
end

@testset "Printing" begin
    C = Component(:C)
    SS = Component(:SS, "Source")
    J0 = Junction(:J0)
    b1 = Bond(C,J0)
    b2 = Bond(J0,SS)
    bg = BondGraph("newbg")
    bgn = BondGraphNode(bg)

    # repr returns the output of the 'show' function
    @test repr(C) == "C:C"
    @test repr(SS) == "SS:Source"
    @test repr(b1) == "Bond C:C ⇀ J0"
    @test repr(b2) == "Bond J0 ⇀ SS:Source"
    @test repr(bg) == "BondGraph newbg (0 Nodes, 0 Bonds)"
    @test repr(bgn) == "BG:newbg"

    add_vertex!(bg, C)
    add_vertex!(bg, SS)
    add_vertex!(bg, J0)
    add_edge!(bg, C, J0)
    add_edge!(bg, J0, SS)
    @test repr(bg) == "BondGraph newbg (3 Nodes, 2 Bonds)"
end

@testset "LightGraph Extra Functions" begin
    c1 = Component(:C)
    c2 = Component(:R)
    c3 = Component(:I)
    j = Junction(:J)

    bg = BondGraph()

    add_vertex!(bg, c1)
    add_vertex!(bg, c2)
    add_vertex!(bg, c3)
    add_vertex!(bg, j)

    add_edge!(bg, c1, j)
    add_edge!(bg, j, c2)
    add_edge!(bg, j, c3)

    # Testing on a selection of common LG functions
    @test Δ(bg) == 3
    @test density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]
end