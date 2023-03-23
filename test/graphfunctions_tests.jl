@testset "BondGraph Properties" begin
    bg = BondGraph("newBG")
    @test name(bg) == "newBG"
    @test isempty(vertices(bg))

    @test eltype(BondGraph) == AbstractNode
    @test eltype(bg) == AbstractNode

    @test isempty(edges(bg))
    @test edgetype(BondGraph) == Graphs.AbstractSimpleEdge{Integer}
    @test edgetype(bg) == Graphs.AbstractSimpleEdge{Integer}
    @test is_directed(bg)

    @test size(zero(BondGraph)) == size(BondGraph())
    @test size(zero(bg)) == size(BondGraph())
end

@testset "Adding and removing elements" begin
    c = Component(:C)
    r = Component(:R)
    j0 = EqualEffort()

    bg = BondGraph()

    @test add_vertex!(bg, c)

    b = add_edge!(bg, (c,1), (j0,1))
    @test src(b) == vertex(c)
    @test dst(b) == vertex(j0)

    add_vertex!(bg, r)
    add_vertex!(bg, j0)

    @test ne(bg) == 1
    @test has_edge(bg, vertex(c), vertex(j0))
    @test !has_edge(bg, vertex(c), vertex(r))

    @test nv(bg) == 3
    @test has_vertex(bg, j0)
    @test has_vertex(bg, 3)
    @test !has_vertex(bg, 0)

    @test components(bg) == [c, r]
    @test junctions(bg) == [j0]

    @test inneighbors(bg, vertex(c)) == []
    @test outneighbors(bg, vertex(c)) == [3]

    @test rem_edge!(bg, c, j0) == b
    @test ne(bg) == 0
    @test rem_vertex!(bg, r)
    @test nv(bg) == 2
end

@testset "BondGraphNode" begin
    bg = BondGraph("RCI")
    bgn = BondGraphNode(bg)

    @test bgn.name == "RCI"
    @test bgn.ports == Bool[]
end

@testset "Printing" begin
    C = Component(:C)
    SS = Component(:SS, "Source")
    J0 = EqualEffort(name="J")
    b1 = Bond(C, J0)
    b2 = Bond(J0, SS)
    port = b1.src
    bg = BondGraph(:newbg)
    bgn = BondGraphNode(bg)

    # repr returns the output of the 'show' function
    @test repr(C) == "C:C"
    @test repr(SS) == "SS:Source"
    @test repr(b1) == "Bond C:C[1] ⇀ J[1]"
    @test repr(b2) == "Bond J[1] ⇀ SS:Source[1]"
    @test repr(port) == "(C:C, 1)"
    @test repr(bg) == "BondGraph newbg (0 Nodes, 0 Bonds)"
    @test repr(bgn) == "BG:newbg"

    add_vertex!(bg, C)
    add_vertex!(bg, SS)
    add_vertex!(bg, J0)
    add_edge!(bg, (C,1), (J0,1))
    add_edge!(bg, (J0,1), (SS,1)) # junction index is '1' here as a quick fix
    @test repr(bg) == "BondGraph newbg (3 Nodes, 2 Bonds)"
end

@testset "Graphs.jl Extra Functions" begin
    c1 = Component(:C)
    c2 = Component(:R)
    c3 = Component(:I)
    j = EqualFlow()

    bg = BondGraph()

    add_vertex!(bg, c1)
    add_vertex!(bg, c2)
    add_vertex!(bg, c3)
    add_vertex!(bg, j)

    add_edge!(bg, (c1,1), (j,1))
    add_edge!(bg, (j,1), (c2,1)) # junction index is '1' here as a quick fix
    add_edge!(bg, (j,1), (c3,1)) # junction index is '1' here as a quick fix

    # Testing on a selection of common graph functions
    @test Δ(bg) == 3
    @test Graphs.density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]
end
