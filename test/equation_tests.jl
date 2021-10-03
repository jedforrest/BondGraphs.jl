@parameters t
D = Differential(t)

@testset "0-junction equations" begin
    model = BondGraph(:RC)
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    connect!(model, R, zero_law)
    connect!(model, zero_law, C)

    @test numports(zero_law) == 2

    @variables E[1:2](t) F[1:2](t)
    @test cr(zero_law) == [
        0 ~ F[1] + F[2],
        0 ~ E[1] - E[2]
    ]
end

@testset "1-junction equations" begin
    c1 = new(:C,:C1)
    c2 = new(:R,:R1)
    c3 = new(:I,:I1)
    j = EqualFlow()

    bg = BondGraph()
    add_node!(bg, [c1, c2, c3, j])
    connect!(bg, c1, j)
    connect!(bg, j, c2)
    connect!(bg, j, c3)

    @test numports(j) == 3
    @test length(j.weights) == 3
    @test j.weights == [1,-1,-1]

    @variables E[1:3](t) F[1:3](t)
    @test cr(j) == [
        0 ~ E[1] - E[2] - E[3],
        0 ~ F[1] + F[2],
        0 ~ F[1] + F[3],
    ]
end

function rlc()
    r = new(:R)
    l = new(:I)
    c = new(:C)
    kvl = EqualEffort(name=:kvl)

    bg = BondGraph()
    add_node!(bg, [c, l, kvl, r])

    connect!(bg, r, kvl)
    connect!(bg, l, kvl)
    connect!(bg, c, kvl)
    return bg
end

@testset "RC circuit" begin
    r = new(:R)
    c = new(:C)
    bg = BondGraph(:RC)
    add_node!(bg, [c, r])
    connect!(bg, r, c)

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)
    @test length(eqs) == 1

    (C,R) = sys.ps
    x = sys.states[1]
    @test eqs == [D(x) ~ -x/C/R]
end

@testset "RLC circuit" begin
    bg = rlc()
    eqs = BondGraphs.equations(bg)
    @test length(eqs) == 2

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)
    (C,L,R) = sys.ps
    (qC,pL) = sys.states
    @test eqs == [
        D(qC) ~ -pL/L - qC/C/R,
        D(pL) ~ qC/C
    ]
end

@testset "Chemical reaction A ⇌ B" begin
    A = new(:ce,:A)
    B = new(:ce,:B)
    re = new(:re,:r)
    bg = BondGraph()

    add_node!(bg,[A,B,re])
    connect!(bg, A, re; dstportindex=1)
    connect!(bg, re, B; srcportindex=2)
    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)

    (xA,xB) = sys.states
    (KA,KB,r) = sys.ps
    @test eqs == [
        D(xA) ~ -r*(KA*xA - KB*xB),
        D(xB) ~ r*(KA*xA - KB*xB)
    ]
end

@testset "Chemical reaction A ⇌ B + C, C ⇌ D" begin
    C_A = new(:ce,:A)
    C_B = new(:ce,:B)
    C_C = new(:ce,:C)
    C_D = new(:ce,:D)
    re1 = new(:re,:r1)
    re2 = new(:re,:r2)
    common_C = EqualEffort()
    BC = EqualFlow()

    bg = BondGraph()
    add_node!(bg,[C_A,C_B,C_C,C_D,re1,re2,common_C,BC])
    connect!(bg,C_A,re1; dstportindex=1)
    connect!(bg,re1,BC; srcportindex=2)
    connect!(bg,BC,C_B)
    connect!(bg,BC,common_C)
    connect!(bg,common_C,C_C)
    connect!(bg,common_C,re2; dstportindex=1)
    connect!(bg,re2,C_D; srcportindex=2)

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)

    (xA,xB,xC,xD) = sys.states
    (KA,KB,KC,KD,r1,r2) = sys.ps
    @test eqs == [
        D(xA) ~ -r1*(KA*xA - KB*xB*KC*xC),
        D(xB) ~ r1*(KA*xA - KB*xB*KC*xC),
        D(xC) ~ r1*(KA*xA - KB*xB*KC*xC) - r2*(KC*xC - KD*xD),
        D(xD) ~ r2*(KC*xC - KD*xD)
    ]
end