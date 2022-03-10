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

@testset "Expose models" begin
    r = Component(:R)
    kcl = EqualFlow(name = :kcl)
    SSA = Component(:SS, :A)
    SSB = Component(:SS, :B)

    bg = BondGraph()
    add_node!(bg, [r, kcl, SSA, SSB])
    connect!(bg, kcl, r)
    connect!(bg, SSA, kcl)
    connect!(bg, kcl, SSB)

    bgn = expose(bg, [SSA, SSB])
    @test numports(bgn) == 2

    eqns = equations(bgn; simplify_eqs = false)
    #TODO: Check whether external variables are in equations
end

@testset "Modular RLC circuit" begin
    r = Component(:R)
    l = Component(:I)
    c = Component(:C)
    kvl = EqualEffort(name = :kvl)
    SS1 = Component(:SS, :SS1)
    SS2 = Component(:SS, :SS2)

    bg1 = BondGraph()
    add_node!(bg1, [r, c, kvl, SS1])
    connect!(bg1, r, kvl)
    connect!(bg1, c, kvl)
    connect!(bg1, SS1, kvl)
    bgn1 = expose(bg1, [SS1])

    bg2 = BondGraph()
    add_node!(bg2, [l, SS2])
    connect!(bg2, l, SS2)
    bgn2 = expose(bg2, [SS2])

    bg = BondGraph()
    add_node(bg, [bgn1, bgn2])
    connect!(bg, bgn1, bgn2)

    eqns = constitutive_relations(bg)
    @test length(eqs) == 2

    sys = ODESystem(bg)
    (C, L, R) = sys.ps
    (qC, pL) = sys.states
    e1 = D(qC) ~ -pL / L + (-qC / C / R)
    e2 = D(pL) ~ qC / C

    @test isequal(simplify(eqs[1].rhs-e1.rhs), 0)
    @test isequal(eqs[2].rhs, e2.rhs)
end

@testset "Modular reaction" begin
    bg1 = BondGraph()
    re = Component(:re, :r)
    SSA = Component(:SS, :SSA)
    SSB = Component(:SS, :SSB)
    add_node!(bg1, [SSA, SSB, re])
    connect!(bg1, SSA, re; dstportindex = 1)
    connect!(bg1, re, SSB; srcportindex = 2)
    bgn1 = expose(bg1, [SSA, SSB])

    bg = BondGraph()
    A = Component(:ce, :A)
    B = Component(:ce, :B)
    add_node!(bg, [A, B, bgn1])
    connect!(bg, A, bgn1; dstportindex = 1)
    connect!(bg, bgn1, B; srcportindex = 2)

    sys = ODESystem(bg)
    eqs = constitutive_relations(bg)

    (xA, xB) = sys.states
    (KA, KB, r) = sys.ps
    e1 = D(xA) ~ -r * (KA * xA - KB * xB)
    e2 = D(xB) ~ r * (KA * xA - KB * xB)

    @test isequal(eqs[1].rhs, e1.rhs)
    @test isequal(eqs[2].rhs, e2.rhs)
end