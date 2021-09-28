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

@testset "1-junction equations" begin
    c1 = new(:C,"C1")
    c2 = new(:R,"R1")
    c3 = new(:I,"I1")
    j = EqualFlow()

    bg = BondGraph()
    add_node!(bg, [c1, c2, c3, j])
    connect!(bg, c1, j)
    connect!(bg, j, c2)
    connect!(bg, j, c3)

    @test j.degree == 3
    @test length(j.weights) == 3
    @test j.weights == [1,-1,-1]

    @variables E_1(t) E_2(t) E_3(t) F_1(t) F_2(t) F_3(t)
    @test BondGraphs.equations(j) == [
        0 ~ E_1 - E_2 - E_3,
        0 ~ F_1 + F_2,
        0 ~ F_1 + F_3,
    ]
end
