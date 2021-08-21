# Based on https://bondgraphtools.readthedocs.io/en/latest/tutorials/RC.html
@testset "BondGraph Construction" begin
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:𝟎)

    add_nodes!(model, [R, C, zero_law])
    @test R in model.nodes
    @test C in model.nodes
    @test zero_law in model.nodes

    connect!(model, R, zero_law)
    connect!(model, zero_law, C)
    @test Bond(R, zero_law) in model.bonds
    @test Bond(zero_law, C) in model.bonds
end

@testset "Construction Failure" begin
    model = BondGraph(name="RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = Junction(:𝟎)

    add_nodes!(model, [R, C, zero_law])
    @test_throws ErrorException add_nodes!(model, R)
    @test_throws ErrorException add_nodes!(model, zero_law)

    connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)
end