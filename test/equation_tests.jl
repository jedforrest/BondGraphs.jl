@parameters t
D = Differential(t)

function RLC()
    r = Component(:R)
    l = Component(:I)
    c = Component(:C)
    kvl = EqualEffort(name = :kvl)

    bg = BondGraph()
    add_node!(bg, [c, l, kvl, r])

    connect!(bg, r, kvl)
    connect!(bg, l, kvl)
    connect!(bg, c, kvl)
    return bg
end

@testset "Equations" begin
    c = Component(:C)
    @parameters C
    @variables E[1](t) F[1](t) q(t)
    @test constitutive_relations(c) == [
        0 ~ q / C - E[1],
        D(q) ~ F[1]
    ]
end

@testset "Parameters" begin
    set_library!(biochemical_library)

    tf = Component(:TF)
    @parameters r
    @test iszero(BondGraphs.params(tf) - [r])

    Ce = Component(:Ce, library = biochemical_library)
    @parameters k R T
    @test iszero(BondGraphs.params(Ce) - [k, R, T])

    Re = Component(:Re, library = biochemical_library)
    @parameters r R T
    @test iszero(BondGraphs.params(Re) - [r, R, T])

    set_library!()
end

@testset "State variables" begin
    r = Component(:R)
    @test isempty(BondGraphs.state_vars(r))

    @variables q(t)
    c = Component(:C)
    @test isequal(BondGraphs.state_vars(c), [q])

    ce = Component(:ce, library = biochemical_library)
    @test isequal(BondGraphs.state_vars(ce), [q])
end

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
    @test constitutive_relations(zero_law) == [
        0 ~ F[1] + F[2],
        0 ~ E[1] - E[2]
    ]
end

@testset "1-junction equations" begin
    c1 = Component(:C, :C1)
    c2 = Component(:R, :R1)
    c3 = Component(:I, :I1)
    j = EqualFlow()

    bg = BondGraph()
    add_node!(bg, [c1, c2, c3, j])
    connect!(bg, c1, j)
    connect!(bg, j, c2)
    connect!(bg, j, c3)

    @test numports(j) == 3
    @test length(j.weights) == 3
    @test j.weights == [1, -1, -1]

    @variables E[1:3](t) F[1:3](t)
    @test constitutive_relations(j) == [
        0 ~ E[1] - E[2] - E[3],
        0 ~ F[1] + F[2],
        0 ~ F[1] + F[3],
    ]
end

@testset "RC circuit" begin
    r = Component(:R)
    c = Component(:C)
    bg = BondGraph(:RC)
    add_node!(bg, [c, r])
    connect!(bg, r, c)

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)
    @test length(eqs) == 1

    (C, R) = sys.ps
    x = sys.states[1]
    e1 = eqs[1]
    e2 = D(x) ~ -x / C / R
    @test e1.lhs == e2.lhs
    @test expand(e1.rhs) == e2.rhs
end

@testset "RLC circuit" begin
    bg = RLC()
    eqs = BondGraphs.equations(bg)
    @test length(eqs) == 2

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)
    (C, L, R) = sys.ps
    (qC, pL) = sys.states
    e1 = D(qC) ~ -pL / L - qC / C / R
    e2 = D(pL) ~ qC / C

    @test expand(eqs[1].rhs) == e1.rhs
    @test eqs[2].rhs == e2.rhs
end

@testset "Chemical reaction A ⇌ B" begin
    set_library!(biochemical_library)

    A = Component(:ce, :A)
    B = Component(:ce, :B)
    re = Component(:re, :r)
    bg = BondGraph()

    add_node!(bg, [A, B, re])
    connect!(bg, A, re; dstportindex = 1)
    connect!(bg, re, B; srcportindex = 2)
    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)

    (xA, xB) = sys.states
    (KA, KB, r) = sys.ps
    e1 = D(xA) ~ -r * (KA * xA - KB * xB)
    e2 = D(xB) ~ r * (KA * xA - KB * xB)

    @test eqs[1].rhs == e1.rhs
    @test eqs[2].rhs == e2.rhs

    set_library!()
end

@testset "Chemical reaction A ⇌ B + C, C ⇌ D" begin
    set_library!(biochemical_library)

    C_A = Component(:ce, :A)
    C_B = Component(:ce, :B)
    C_C = Component(:ce, :C)
    C_D = Component(:ce, :D)
    re1 = Component(:re, :r1)
    re2 = Component(:re, :r2)
    common_C = EqualEffort()
    BC = EqualFlow()

    bg = BondGraph()
    add_node!(bg, [C_A, C_B, C_C, C_D, re1, re2, common_C, BC])
    connect!(bg, C_A, re1; dstportindex = 1)
    connect!(bg, re1, BC; srcportindex = 2)
    connect!(bg, BC, C_B)
    connect!(bg, BC, common_C)
    connect!(bg, common_C, C_C)
    connect!(bg, common_C, re2; dstportindex = 1)
    connect!(bg, re2, C_D; srcportindex = 2)

    sys = ODESystem(bg)
    eqs = ModelingToolkit.equations(sys)

    (xA, xB, xC, xD) = sys.states
    (KA, KB, KC, KD, r1, r2) = sys.ps
    e1 = D(xA) ~ -r1 * (KA * xA - KB * xB * KC * xC)
    e2 = D(xB) ~ r1 * (KA * xA - KB * xB * KC * xC)
    e3 = D(xC) ~ r1 * (KA * xA - KB * xB * KC * xC) - r2 * (KC * xC - KD * xD)
    e4 = D(xD) ~ r2 * (KC * xC - KD * xD)

    # equations are not simplifying with exp/log rules
    @test_broken eqs[1].rhs == e1.rhs
    @test_broken eqs[2].rhs == e2.rhs
    @test_broken eqs[3].rhs == e3.rhs
    @test eqs[4].rhs == e4.rhs

    set_library!()
end