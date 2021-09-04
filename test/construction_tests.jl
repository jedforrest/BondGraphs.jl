# Based on https://bondgraphtools.readthedocs.io/en/latest/tutorials/RC.html
@testset "BondGraph Construction" begin
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:J0)

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
    model = BondGraph(name="RCI")
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = Junction(:J0)
    one_law = Junction(:J1)

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
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:J0)

    add_node!(model, [R, C, zero_law])
    @test_throws ErrorException add_node!(model, R)
    @test_throws ErrorException add_node!(model, zero_law)

    connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)

    one_law = Junction(:J1)
    @test_throws ErrorException remove_node!(model, one_law)
    @test_throws ErrorException swap!(model, C, one_law)
end