@parameters t
D = Differential(t)

function RLC()
    r = Component(:R)
    l = Component(:I)
    c = Component(:C)
    kvl = EqualEffort(name=:kvl)

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
    cr = [
        0 ~ q / C - E[1],
        D(q) ~ F[1]
    ]

    @test isequal(equations(c), cr)
    @test isequal(constitutive_relations(c), cr)

    j = EqualEffort()
    @test isequal(equations(j), Equation[])
end

@testset "Parameters" begin
    tf = Component(:TF)
    @parameters n
    @test iszero(parameters(tf) - [n])

    Ce = Component(:Ce)
    @parameters K
    @test iszero(parameters(Ce) - [K])

    bg = RLC()
    @parameters C L R
    @test iszero(parameters(bg) - [C, L, R])
end

@testset "Globals" begin
    re = Component(:Re)
    c = Component(:C)
    @parameters R T

    @test iszero(globals(re) - [T, R])
    @test globals(c) == Num[]
end

@testset "State variables" begin
    r = Component(:R)
    @test isempty(states(r))

    @variables q(t)
    c = Component(:C)
    @test isequal(states(c), [q])

    ce = Component(:ce)
    @test isequal(states(ce), [q])

    bg = RLC()
    @variables q(t) p(t)
    @test iszero(states(bg) - [q, p])
end

@testset "Controls" begin
    se = Component(:Se)
    sf = Component(:Sf)
    c = Component(:C)
    @parameters fs(t) es(t)

    @test isequal(controls(se), [es])
    @test isequal(controls(sf), [fs])
    @test isequal(controls(c), Num[])
end

@testset "Constitutive relations" begin
    eqE = EqualEffort()
    eqF = EqualFlow()
    @test constitutive_relations(eqE) == Equation[]
    @test constitutive_relations(eqF) == Equation[]

    bg = RLC()
    cr_bg = constitutive_relations(bg)
    sys = ODESystem(bg)

    C, L, R = sys.ps
    q, p = sys.states
    cr1 = D(q) ~ -(q / C) / R + (-p) / L
    cr2 = D(p) ~ q / C

    @test isequal(cr_bg[1].lhs, cr1.lhs)
    @test_broken isequal(simplify(cr_bg[1].rhs - cr1.rhs), 0)
    @test isequal(cr_bg[2], cr2)

    cr_bgn = constitutive_relations(BondGraphNode(bg))
    @test isequal(cr_bgn[1].lhs, cr1.lhs)
    @test_broken isequal(simplify(cr_bgn[1].rhs - cr1.rhs), 0)
    @test isequal(cr_bgn[2], cr2)
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
    @test isequal(constitutive_relations(zero_law), [
        0 ~ F[1] + F[2],
        0 ~ E[1] - E[2]
    ])
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
    @test isequal(constitutive_relations(j), [
        0 ~ E[1] - E[2] - E[3],
        0 ~ F[1] + F[2],
        0 ~ F[1] + F[3],
    ])
end

@testset "RC circuit" begin
    r = Component(:R)
    c = Component(:C)
    bg = BondGraph(:RC)
    add_node!(bg, [c, r])
    connect!(bg, r, c)

    sys = ODESystem(bg)
    eqs = constitutive_relations(bg)
    @test_broken length(eqs) == 1

    (C, R) = sys.ps
    x = sys.states[1]
    e1 = eqs[1]
    e2 = D(x) ~ -x / C / R

    @test isequal(e1.lhs, e2.lhs)
    @test_broken isequal(expand(e1.rhs), e2.rhs)
end

@testset "RL circuit" begin
    r = Component(:R)
    l = Component(:I)
    bg = BondGraph(:RL)
    add_node!(bg, [r, l])
    connect!(bg, l, r)

    eqs = constitutive_relations(bg)
    sys = ODESystem(bg)
    x = sys.states[1]
    (R, L) = sys.ps
    @test eqs == [D(x) ~ -R * x / L]
end

@testset "RLC circuit" begin
    bg = RLC()
    eqs = constitutive_relations(bg)
    @test_broken length(eqs) == 2

    sys = ODESystem(bg)
    (C, L, R) = sys.ps
    (qC, pL) = sys.states
    e1 = D(qC) ~ -pL / L + (-qC / C / R)
    e2 = D(pL) ~ qC / C

    @test_broken isequal(simplify(eqs[1].rhs - e1.rhs), 0)
    @test isequal(eqs[2].rhs, e2.rhs)
end

@testset "Chemical reaction A ⇌ B" begin
    A = Component(:ce, :A)
    B = Component(:ce, :B)
    re = Component(:re, :r)
    bg = BondGraph()

    add_node!(bg, [A, B, re])
    connect!(bg, A, re; dstportindex=1)
    connect!(bg, re, B; srcportindex=2)
    sys = ODESystem(bg)
    eqs = constitutive_relations(bg)

    (xA, xB) = sys.states
    (KA, KB, r) = sys.ps
    e1 = D(xA) ~ -r * (KA * xA - KB * xB)
    e2 = D(xB) ~ r * (KA * xA - KB * xB)

    @test isequal(eqs[1].rhs, e1.rhs)
    @test isequal(eqs[2].rhs, e2.rhs)
end

@testset "Chemical reaction A ⇌ B + C, C ⇌ D" begin
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
    connect!(bg, C_A, re1; dstportindex=1)
    connect!(bg, re1, BC; srcportindex=2)
    connect!(bg, BC, C_B)
    connect!(bg, BC, common_C)
    connect!(bg, common_C, C_C)
    connect!(bg, common_C, re2; dstportindex=1)
    connect!(bg, re2, C_D; srcportindex=2)

    sys = ODESystem(bg)
    eqs = constitutive_relations(bg)

    (xA, xB, xC, xD) = sys.states
    (KA, KB, KC, KD, r1, r2) = sys.ps
    e1 = D(xA) ~ -r1 * (KA * xA - KB * xB * KC * xC)
    e2 = D(xB) ~ r1 * (KA * xA - KB * xB * KC * xC)
    e3 = D(xC) ~ r1 * (KA * xA - KB * xB * KC * xC) - r2 * (KC * xC - KD * xD)
    e4 = D(xD) ~ r2 * (KC * xC - KD * xD)

    # equations are not simplifying with exp/log rules
    @test isequal(eqs[1].rhs, e1.rhs)
    @test isequal(eqs[2].rhs, e2.rhs)
    @test isequal(eqs[3].rhs, e3.rhs)
    @test isequal(eqs[4].rhs, e4.rhs)
end

