@testset "SS component system" begin
    SS = SourceSensor(name = :SS)

    @test length(freeports(SS)) == 1
    @test numports(SS) == 1
    @test length(parameters(SS)) == 0
    @test length(states(SS)) == 0
    @test length(equations(SS)) == 0
    @test length(constitutive_relations(SS)) == 0

    # Todo: fix tests
    sys = ODESystem(SS)
    @test length(sys.systems) == 1
    @test sys.p1.E isa Num
    @test sys.p1.F isa Num
end

function find_subsys(sys,s)
    subsys = ModelingToolkit.get_systems(sys)
    return filter(x -> nameof(x) == s, subsys)[1]
end

@testset "Expose models" begin
    r = Component(:R)
    kcl = EqualFlow(name = :kcl)
    SSA = SourceSensor(name = :A)
    SSB = SourceSensor(name = :B)

    bg = BondGraph()
    add_node!(bg, [r, kcl, SSA, SSB])
    connect!(bg, kcl, r)
    connect!(bg, SSA, kcl)
    connect!(bg, kcl, SSB)

    bgn = expose(bg, [SSA, SSB])
    @test numports(bgn) == 2

    sys = ODESystem(bgn) 
    expanded_sys = expand_connections(sys) # Note that the equations shouldn't simplify by default
    eqns = equations(expanded_sys)
    
    p1 = find_subsys(sys,:p1)
    p2 = find_subsys(sys,:p2)
    (E1,F1,E2,F2) = (p1.E, p1.F, p2.E, p2.F)

    Asys = find_subsys(sys,:A)
    Bsys = find_subsys(sys,:B)
    (AE,AF,BE,BF) = (Asys.p1.E, Asys.p1.F, Bsys.p1.E, Bsys.p1.F)

    @test (E1 ~ AE) in eqns
    @test (F1 ~ AF) in eqns
    @test (E2 ~ BE) in eqns
    @test (F2 ~ BF) in eqns
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