@parameters t
D = Differential(t)

@testset "0-junction equations" begin
    model = BondGraph("RC")
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    connect!(model, R, zero_law)
    connect!(model, zero_law, C)

    @test zero_law.degree == 2

    @variables E_1(t) E_2(t) F_1(t) F_2(t)
    @test BondGraphs.equations(zero_law) == [
        0 ~ F_1 + F_2,
        0 ~ E_1 - E_2
    ]
end