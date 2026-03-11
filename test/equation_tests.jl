@variables e(t) f(t) q(t)
@parameters C L R

@testset "Component Systems" begin
    eqs = [D(q) ~ f, q ~ C * e]
    @named cap = StaticStorageElement(eqs, [e], [f], [q])

    sys = system(cap)
    @test nameof(sys) == :cap
    @test length(equations(cap)) == 4  # includes port relations
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
    @test isempty(equations(bg))

    @named cap = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q])
    add_comp!(bg, cap)
    @test equations(bg) == [D(cap.sys.q) ~ 0.0]
end

@testset "Bond Graph RC System" begin
    # System 1
    using BondGraphs: capacitor, resistor, KVL
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
    @test isequal(cr_bg[2], cr2)

    # With default values
    cr_subbed = constitutive_relations(bg; sub_defaults=true)
    @test isequal(cr_subbed, [D(q) ~ p, D(p) ~ -1 - p - q])
end

# TODO CONTINUE FROM HERE
# @testset "Chemical reaction A ⇌ B" begin
#     A = Component(:ce, :A)
#     B = Component(:ce, :B)
#     re = Component(:re, :r)
#     bg = BondGraph()

#     add_node!(bg, [A, B, re])
#     connect!(bg, A, (re, 1))
#     connect!(bg, (re, 2), B)
#     sys = ODESystem(bg)
#     eqs = sorted_eqs(sys)

#     (xA, xB) = (sys.A.q, sys.B.q)
#     (KA, KB, r) = (sys.A.K, sys.B.K, sys.r.r)
#     e1 = D(xA) ~ r * (-KA * xA + KB * xB)
#     e2 = D(xB) ~ r * (KA * xA - KB * xB)

#     @test isequal(eqs[1].rhs, e1.rhs)
#     @test isequal(eqs[2].rhs, e2.rhs)
# end

# @testset "Chemical reaction A ⇌ B + C, C ⇌ D" begin
#     C_A = Component(:ce, :A)
#     C_B = Component(:ce, :B)
#     C_C = Component(:ce, :C)
#     C_D = Component(:ce, :D)
#     re1 = Component(:re, :r1)
#     re2 = Component(:re, :r2)
#     common_C = EqualEffort()
#     BC = EqualFlow()

#     bg = BondGraph()
#     add_node!(bg, [C_A, C_B, C_C, C_D, re1, re2, common_C, BC])
#     connect!(bg, C_A, (re1, 1))
#     connect!(bg, (re1, 2), BC)
#     connect!(bg, BC, C_B)
#     connect!(bg, BC, common_C)
#     connect!(bg, common_C, C_C)
#     connect!(bg, common_C, (re2, 1))
#     connect!(bg, (re2, 2), C_D)

#     sys = ODESystem(bg)
#     eqs = sorted_eqs(sys)

#     (xA, xB, xC, xD) = (sys.A.q, sys.B.q, sys.C.q, sys.D.q)
#     (KA, KB, KC, KD, r1, r2) = (sys.A.K, sys.B.K, sys.C.K, sys.D.K, sys.r1.r, sys.r2.r)
#     e1 = D(xA) ~ -r1 * (KA * xA - KB * xB * KC * xC)
#     e2 = D(xB) ~ r1 * (KA * xA - KB * xB * KC * xC)
#     e3 = D(xC) ~ r1 * (KA * xA - KB * xB * KC * xC) - r2 * (KC * xC - KD * xD)
#     e4 = D(xD) ~ r2 * (KC * xC - KD * xD)

#     @test isequal(simplify(eqs[1].rhs - e1.rhs), 0)
#     @test isequal(simplify(eqs[2].rhs - e2.rhs), 0)
#     @test isequal(simplify(eqs[3].rhs - e3.rhs), 0)
#     @test isequal(simplify(eqs[4].rhs - e4.rhs), 0)
# end
