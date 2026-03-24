@variables e(t) f(t) q(t)
@parameters C L R

@testset "Component Systems" begin
    eqs = [D(q) ~ f, q ~ C * e]
    @named cap = StaticStorageElement(eqs, [e], [f], [q])

    sys = system(cap)
    @test nameof(sys) == :cap
    @test length(equations(expand_connections(cap.sys))) == 3  # includes port relations
    @test length(constitutive_relations(cap)) == 2
    @test isequal(constitutive_relations(cap), eqs)

    j = EqualEffort(; name=:𝟎)
    @test nameof(system(j)) == :𝟎
    @test isempty(equations(j))
    @test length(constitutive_relations(j)) == 1

    tf = transformer(2)
    @test nameof(system(tf)) == Symbol("TF{2}")
    @test length(equations(tf)) == 2
end

@testset "Bond Graph System" begin
    bg = BondGraph()

    @test isnothing(bg.sys)
    compile_system!(bg)
    @test isempty(equations(bg.sys))

    @named cap = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q])
    @named res = DissipatorElement([e ~ R * f], [e], [f])
    add_comp!(bg, cap)
    add_comp!(bg, res)
    connect!(bg, cap, res)
    @test full_equations(system(bg)) == [D(cap.sys.q) ~ -cap.sys.q / (cap.sys.C*res.sys.R)]
end

@testset "Bond Graph RC System" begin
    # System 1
    @named cap = capacitor()
    @named res = resistor()
    @named kvl = KVL()
    bg = BondGraph([cap, res, kvl])
    connect!(bg, cap, kvl)
    connect!(bg, res, kvl)

    sys = system(bg, simplify=false)
    @test length(full_equations(sys)) == 11

    sys = system(bg)
    @test length(full_equations(sys)) == 1

    cr = constitutive_relations(bg)
    q, C, R = sys.cap.q, sys.cap.C, sys.res.R
    @test isequal(cr, [D(q) ~ -q / (C * R)])
end

@testset "Bond Graph RCI circuit" begin
    # System 2
    bg = RCI()
    cr_bg = constitutive_relations(bg)
    sys = system(bg)

    C, L, R, E = (sys.c.C, sys.i.L, sys.r.R, sys.v.E)
    q, p = sys.c.q, sys.i.λ

    cr1 = D(q) ~ p / L
    cr2 = D(p) ~ -E + (-q) / C + (-R*p) / L

    # Constitutive relations
    @test isequal(cr_bg[1], cr1)
    @test isequal(cr_bg[2], simplify(cr2))

    # With default values
    cr_subbed = constitutive_relations(bg; sub_defaults=true)
    @test isequal(cr_subbed, [D(q) ~ p, D(p) ~ -1 - p - q])
end

@testset "Equation Rewriter" begin
    # simplifying exp/log for chemical bond graphs
    rewriter = BondGraphs.rewriter
    @variables a b

    expr = 2log(a) + log(b)
    simplified_expr = simplify(expr; rewriter)
    @test isequal(simplified_expr, log((b * a^2)))

    expr = -2b*exp(4log(a)) + 3
    simplified_expr = simplify(expr; rewriter)
    @test isequal(simplified_expr, 3 - 2b*a^4)

    @variables k R T K1 K2 K3 K4 x1 x2 x3 x4
    expr = k*(exp(log(K4*x4)) - exp((R*T*log(K1*x1) + R*T*log(K3*x3) + R*T*log(K2*x2)) / (R*T)))
    simplified_expr = simplify(simplify(expr); rewriter)  # double simplification needed
    @test isequal(simplified_expr, k * (K4*x4 - K1*K2*K3*x1*x2*x3))
end

@testset "Chemical reaction A ⇌ B" begin
    @named A = chemicalspecies()
    @named B = chemicalspecies()
    @named re = reaction()
    bg = BondGraph()

    for c in [A, B, re]
        add_comp!(bg, c)
    end
    connect!(bg, A, re[1])
    connect!(bg, re[2], B)

    sys = system(bg)
    xA, xB = unknowns(sys)
    KB, KA, k = parameters(sys)

    xA, xB = sys.A.x, sys.B.x
    KA, KB, k = sys.A.K, sys.B.K, sys.re.κ

    cr = constitutive_relations(bg)
    @test isequal(cr[1], D(xA) ~ k*(-KA*xA + KB*xB))
    @test isequal(cr[2], D(xB) ~ k*(KA*xA - KB*xB))
end

@testset "EqualFlow equations" begin
    bg = begin
        @named J = EqualFlow()
        @named C1 = chemicalspecies()
        @named C2 = chemicalspecies()
        @named C3 = chemicalspecies()
        @named C4 = chemicalspecies()
        @named Re = reaction()
        BondGraph([J, C1, C2, C3, C4, Re])
    end

    # no connections
    @test BondGraphs.port_weight.(J.ports) == [0]

    # C1, C2, C3 -> J -> Re -> C4
    connect!(bg, C1, J)
    connect!(bg, C2, J)
    connect!(bg, C3, J)
    connect!(bg, J, Re)
    connect!(bg, Re, C4)
    @test BondGraphs.port_weight.(J.ports) == [1, 1, 1, -1]
    @test BondGraphs.port_weight.(Re.ports) == [1, -1]
    cr1 = constitutive_relations(bg)

    # C1, C2, C3 <- J -> Re -> C4
    disconnect!(bg, C1, J)
    disconnect!(bg, C2, J)
    disconnect!(bg, C3, J)
    connect!(bg, J, C1)
    connect!(bg, J, C2)
    connect!(bg, J, C3)
    @test BondGraphs.port_weight.(J.ports) == [-1, -1, -1, -1]
    cr2 = constitutive_relations(bg)

    # check that the CR are the same, regardless of internal bond directions
    @test all(isequal.(cr1, cr2))
end

# TODO CONTINUE FROM HERE
# Fix simplification rules
@testset "Chemical reaction A ⇌ B + C, C ⇌ E" begin
    @named A = chemicalspecies()
    @named B = chemicalspecies()
    @named C = chemicalspecies()
    @named E = chemicalspecies()
    @named re1 = reaction()
    @named re2 = reaction()
    @named common_C = EqualEffort()
    @named BC = EqualFlow()
    bg = BondGraph()

    for c in [A, B, C, E, re1, re2, common_C, BC]
        add_comp!(bg, c)
    end
    connect!(bg, A, re1[1])
    connect!(bg, re1[2], BC)
    connect!(bg, BC, B)
    connect!(bg, BC, common_C)
    connect!(bg, common_C, C)
    connect!(bg, common_C, re2[1])
    connect!(bg, re2[2], E)

    sys = system(bg)

    observed(sys)

    constitutive_relations(bg[:BC])

    # TODO CONTINUE FROM HERE
    BondGraphs.inneighbor_comps(bg, bg[:BC])
    BondGraphs.outneighbor_comps(bg, bg[:BC])
    bg[:BC].ports

    full_equations(sys)
    cr = constitutive_relations(bg)

    (xA, xB, xC, xE) = (sys.A.x, sys.B.x, sys.C.x, sys.E.x)
    (KA, KB, KC, KE, k1, k2) = (sys.A.K, sys.B.K, sys.C.K, sys.E.K, sys.re1.κ, sys.re2.κ)
    eq1 = D(xA) ~ -k1 * (KA * xA - KB * xB * KC * xC)
    eq2 = D(xB) ~ k1 * (KA * xA - KB * xB * KC * xC)
    eq3 = D(xC) ~ k1 * (KA * xA - KB * xB * KC * xC) - k2 * (KC * xC - KE * xE)
    eq4 = D(xE) ~ k2 * (KC * xC - KE * xE)

    @test isequal(cr[1], eq1)
    @test isequal(cr[2], eq2)
    @test isequal(cr[3], eq3)
    @test isequal(cr[4], eq4)
end
