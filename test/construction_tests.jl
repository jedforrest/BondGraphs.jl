function RCI()
    model = BondGraph("RCI")
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = Junction(:𝟎)

    add_node!(model, [C, R, I, SS, zero_law])
    connect!(model, R, zero_law)
    connect!(model, C, zero_law)
    connect!(model, zero_law, I)
    connect!(model, zero_law, SS)

    model
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

    disconnect!(model, zero_law, I)
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
    C = new(:C)
    R = new(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    @test_logs (:warn, "R:R already in model") add_node!(model, R)
    @test_logs (:warn, "J0 already in model") add_node!(model, zero_law)

    connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)

    tf = new(:TF)
    @test_throws ErrorException remove_node!(model, tf)
    @test_throws ErrorException swap!(model, C, tf)
    
    one_law = Junction(:J1)
    @test_logs (:warn, "J1 not in model") remove_node!(model, one_law)
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
    tf = new(:TF,:n)
    @test tf isa Component{2}
    @test tf.type == :TF
    @test numports(tf) == 2

    r = new(:R)
    @test r isa Component{1}
    @test r.type == :R
end

@parameters t
D = Differential(t)

@testset "Equations" begin
    c = new(:C)
    @parameters C
    @variables E[1](t) F[1](t) q(t)
    @test cr(c) == [
        0 ~ q/C - E[1],
        D(q) ~ F[1]
    ]
end

@testset "Parameters" begin
    tf = new(:TF)
    @parameters r
    @test iszero(params(tf) - [r])

    Ce = new(:Ce)
    @parameters k R T
    @test iszero(params(Ce) - [k, R, T])

    Re = new(:Re)
    @parameters r R T
    @test iszero(params(Re) - [r, R, T])
end

@testset "State variables" begin
    r = new(:R)
    @test isempty(state_vars(r))

    @variables q(t)
    c = new(:C)
    @test isequal(state_vars(c), [q])

    ce = new(:ce)
    @test isequal(state_vars(c), [q])
end
@testset "Inserting Nodes" begin
    bg = RCI()

    c, r, 𝟎 = bg.nodes[[1,2,5]]

    bondc0 = getbonds(bg, c, 𝟎)[1]
    bondr0 = getbonds(bg, r, 𝟎)[1]

    tf = Component(:TF, numports=2)
    insert_node!(bg, bondc0, tf)
    insert_node!(bg, bondr0, Junction(:𝟏))

    @test tf in bg.nodes
    @test nv(bg) == 7
    @test ne(bg) == 6
end

@testset "Merging components" begin
    bg = RCI()

    newC = Component(:C, "newC")
    newR = Component(:R, "newR")
    add_node!(bg, [newC, newR])
    connect!(bg, newC, newR)

    C = getnodes(bg, "C")[1]
    merge_nodes!(bg, C, newC)

    R = getnodes(bg, "R")[1]
    merge_nodes!(bg, R, newR; junction=Junction(:𝟏))

    @test isempty(getnodes(bg, "newC"))
    @test length(getnodes(bg, :𝟏)) == 1
    @test length(getnodes(bg, :𝟎)) == 2
    @test nv(bg) == 7
    @test ne(bg) == 7
end

@testset "Simplifying Junctions" begin
    bg = RCI()
    C, R, I, SS, 𝟎 = bg.nodes

    J0_new_1 = Junction(:𝟎, "new0_1")
    J0_new_2 = Junction(:𝟎, "new0_2")
    insert_node!(bg, (C, 𝟎), J0_new_1)
    insert_node!(bg, (R, 𝟎), J0_new_2)
    connect!(bg, J0_new_1, J0_new_2)

    J1_new_1 = Junction(:𝟏)
    J1_new_2 = Junction(:𝟏)
    add_node!(bg, J1_new_1)
    connect!(bg, 𝟎, J1_new_1)
    insert_node!(bg, (SS, 𝟎), J1_new_2)

    # Removing junction redundancies
    simplify_junctions!(bg, squash_identical=false)
    @test length(getnodes(bg, :𝟏)) == 0
    @test nv(bg) == 7
    @test ne(bg) == 7

    # Squashing junction duplicates into a single junction
    simplify_junctions!(bg)
    @test length(getnodes(bg, :𝟎)) == 1
    @test nv(bg) == 5
    @test ne(bg) == 4
end
