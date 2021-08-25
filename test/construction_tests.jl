# Based on https://bondgraphtools.readthedocs.io/en/latest/tutorials/RC.html
@testset "BondGraph Construction" begin
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:𝟎)

    add_node!(model, [R, C, zero_law])
    @test R in model.nodes
    @test C in model.nodes
    @test zero_law in model.nodes

    connect!(model, R, zero_law)
    connect!(model, zero_law, C)
    @test Bond(R, zero_law) in model.bonds
    @test Bond(zero_law, C) in model.bonds
end

@testset "BondGraph Modification" begin
    model = BondGraph(name="RCI")
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = Junction(:𝟎)
    one_law = Junction(:𝟏)

    add_node!(model, [C, R, I, SS, zero_law, one_law])
    remove_node!(model, [SS, one_law])
    @test !(SS in model.nodes)
    @test !(one_law in model.nodes)

    connect!(model, R, zero_law)
    connect!(model, C, zero_law)
    connect!(model, zero_law, I)
    @test ne(model) == 3
    @test Bond(zero_law, I) in model.bonds

    disconnect!(model, zero_law, I)
    @test ne(model) == 2
    @test !(Bond(zero_law, I) in model.bonds)

    connect!(model, zero_law, I)
    swap!(model, zero_law, one_law)
    @test one_law in model.nodes
    @test inneighbors(model, one_law) == [R, C]
    @test outneighbors(model, one_law) == [I]
end

@testset "Construction Failure" begin
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:𝟎)

    add_node!(model, [R, C, zero_law])
    @test_throws ErrorException add_node!(model, R)
    @test_throws ErrorException add_node!(model, zero_law)

    connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)

    one_law = Junction(:𝟏)
    @test_throws ErrorException remove_node!(model, one_law)
    @test_throws ErrorException swap!(model, C, one_law)
end