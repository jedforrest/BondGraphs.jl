function RCI(name=:RCI)
    model = BondGraph(name)
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = EqualEffort()

    add_node!(model, [C, R, I, SS, zero_law])
    connect!(model, R, zero_law)
    connect!(model, C, zero_law)
    connect!(model, zero_law, I)
    connect!(model, zero_law, SS)

    model
end

@testset "Creating Components" begin
    C = Component(:C)
    @test type(C) == "C"
    @test name(C) == "C"

    R = Component(:R, "newR")
    @test name(R) == "newR"

    I = Component(:I; L=5)
    @test I.L == 5

    C2 = Component(:C; q=2)
    @test C2.q == 2

    SS = SourceSensor()
    @test type(SS) == "SS"
end

@testset "Creating Junctions" begin
    EqE_1 = EqualEffort()
    EqE_2 = EqualEffort(name="foo")
    EqF = EqualFlow()

    @test name(EqE_1) == "𝟎"
    @test name(EqE_2) == "foo"
    @test name(EqF) == "𝟏"
end

# Based on https://bondgraphtools.readthedocs.io/en/latest/tutorials/RC.html
@testset "BondGraph Construction" begin
    model = BondGraph(:RC)
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    @test R in model.nodes
    @test C in model.nodes
    @test zero_law in model.nodes

    b1 = connect!(model, R, zero_law)
    b2 = connect!(model, zero_law, C)
    @test b1 in model.bonds
    @test b2 in model.bonds
end

@testset "BondGraph Modification" begin
    model = BondGraph(:RCI)
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = EqualEffort()
    one_law = EqualFlow()

    add_node!(model, [C, R, I, SS, zero_law, one_law])
    remove_node!(model, [SS, one_law])
    @test !(SS in model.nodes)
    @test !(one_law in model.nodes)

    connect!(model, R, zero_law)
    connect!(model, C, zero_law)

    @test I.freeports == [true]
    b1 = connect!(model, zero_law, I)
    @test ne(model) == 3
    @test b1 in model.bonds
    @test I.freeports == [false]

    disconnect!(model, I, zero_law) # test disconnect when node order is swapped
    @test ne(model) == 2
    @test !(b1 in model.bonds)
    @test I.freeports == [true]

    connect!(model, zero_law, I)
    swap!(model, zero_law, one_law)
    @test I.freeports == [false]
    @test one_law in model.nodes
    @test inneighbors(model, one_law) == [R, C]
    @test outneighbors(model, one_law) == [I]
end

@testset "Construction Failure" begin
    model = BondGraph(:RC)
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    @test_logs (:warn, "Node 'R' already in model") add_node!(model, R)
    @test_logs (:warn, "Node '𝟎_3' already in model") add_node!(model, zero_law)

    bond = connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)

    one_law = EqualFlow()
    @test_logs (:warn, "Node '𝟏' not in model") remove_node!(model, one_law)

    tf = Component(:TF)
    add_node!(model, tf)
    @test_throws ErrorException swap!(model, tf, C)

    # if inserting a node fails, the original nodes should still remain connected
    @test_throws ErrorException insert_node!(model, bond, Component(:I))
    @test has_edge(model, bond)
end

@testset "Chemical reaction" begin
    model = BondGraph(:Chemical)
    A = Component(:C, :A)
    B = Component(:C, :B)
    C = Component(:C, :C)
    D = Component(:C, :D)
    Re = Component(:Re, :Reaction, numports=2)
    J_AB = EqualFlow()
    J_CD = EqualFlow()

    add_node!(model, [A, B, C, D, Re, J_AB, J_CD])
    connect!(model, A, J_AB)
    connect!(model, B, J_AB)
    connect!(model, C, J_CD)
    connect!(model, D, J_CD)

    @test freeports(Re) == [true, true]
    @test freeports(J_AB) == [false, false]

    # Connecting junctions to specific ports in Re
    connect!(model, Re, J_CD, srcportindex=2)
    @test freeports(Re) == [true, false]

    # connecting to a full port should fail
    @test_throws ErrorException connect!(model, J_AB, Re, dstportindex=2)

    connect!(model, J_AB, Re, dstportindex=1)
    @test freeports(Re) == [false, false]

    @test nv(model) == 7
    @test ne(model) == 6
end

@testset "Standard components" begin
    tf = Component(:TF, :n)
    @test tf isa Component{2}
    @test tf.type == "TF"
    @test numports(tf) == 2

    r = Component(:R)
    @test r isa Component{1}
    @test r.type == "R"
end

@testset "Inserting Nodes" begin
    bg = RCI()

    c, r, J0 = bg.nodes[[1, 2, 5]]

    bondc0 = getbonds(bg, c, J0)[1]
    bondr0 = getbonds(bg, r, J0)[1]

    tf = Component(:TF, numports=2)
    insert_node!(bg, bondc0, tf)
    insert_node!(bg, bondr0, EqualFlow())

    @test tf in bg.nodes
    @test nv(bg) == 7
    @test ne(bg) == 6
end

@testset "Merging components" begin
    bg = RCI()
    C = bg.C
    R = bg.R

    newC = Component(:C, :newC)
    newR = Component(:R, :newR)
    add_node!(bg, [newC, newR])
    connect!(bg, newC, newR)

    merge_nodes!(bg, C, newC)
    merge_nodes!(bg, R, newR; junction=EqualFlow())

    @test isempty(getnodes(bg, "newC"))
    @test length(getnodes(bg, EqualFlow)) == 1
    @test length(getnodes(bg, EqualEffort)) == 2
    @test nv(bg) == 7
    @test ne(bg) == 7
end

@testset "Simplifying Junctions" begin
    bg = RCI()
    C, R, I, SS, J0 = bg.nodes

    J0_new_1 = EqualEffort(; name=:new0_1)
    J0_new_2 = EqualEffort(; name=:new0_2)
    insert_node!(bg, (C, J0), J0_new_1)
    insert_node!(bg, (R, J0), J0_new_2)
    connect!(bg, J0_new_1, J0_new_2)

    J1_new_1 = EqualFlow(; name=:new1_1)
    J1_new_2 = EqualFlow(; name=:new1_2)
    add_node!(bg, J1_new_1)
    connect!(bg, J0, J1_new_1)
    insert_node!(bg, (SS, J0), J1_new_2)

    # Removing junction redundancies
    @test length(getnodes(bg, EqualFlow)) == 2
    simplify_junctions!(bg, squash_identical=false)
    @test length(getnodes(bg, EqualFlow)) == 0
    @test nv(bg) == 7
    @test ne(bg) == 7

    # Squashing junction duplicates into a single junction
    simplify_junctions!(bg)
    @test length(getnodes(bg, EqualEffort)) == 1
    @test nv(bg) == 5
    @test ne(bg) == 4
end

@testset "BondGraphNodes" begin
    C = Component(:C, "C")
    bg1 = BondGraph("first")
    bg2 = BondGraph("second")
    bg3 = BondGraph("third")
    main = BondGraph("Main")

    bgn1 = BondGraphNode(bg1)
    bgn2 = BondGraphNode(bg2)
    bgn3 = BondGraphNode(bg3)

    @test bgn1.bondgraph === bg1
    @test bgn1.type === "BG"
    @test bgn1.name === bg1.name
    @test freeports(bgn1) == Bool[]

    add_node!(bg1, C)
    add_node!(bg2, bgn1)
    add_node!(bg3, bgn2)
    add_node!(main, bgn3)

    @test main.third.second.first.C === C

    C2 = Component(:C, "C") # Same name
    add_node!(bg1, C2)
    @test main.third.second.first.C == [C, C2]
end
