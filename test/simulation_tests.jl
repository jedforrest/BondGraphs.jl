@testset "Setting variables" begin
    c = Component(:C)
    re = Component(:Re)

    # getting
    @test c.C == 1 && c.q == 0
    @test re.r == 1 && re.R ≈ 8.314 && re.T == 310

    # setting
    c.C = 2
    re.T = 200
    @test c.C == 2
    @test re.T == 200
end

@testset "Setting non-numeric control variables" begin
    f(t) = sin(2t) # forcing function

    sf = Component(:Sf)
    @test sf.fs(1) ≈ 1

    sf.fs = f
    @test sf.fs == f
    @test sf.fs(1) ≈ f(1)
end

@testset "Incompatible variables fail" begin
    re = Component(:Re)
    c = Component(:C)
    @test_throws ErrorException re.s = 1
    @test_throws ErrorException c.p = 2
end

@testset "Simulate RC circuit" begin
    r = Component(:R; R=2)
    c = Component(:C; C=1, q=10)
    bg = BondGraph(:RC)

    add_node!(bg, [c, r])
    connect!(bg, r, c)

    f(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sol = simulate(bg, tspan)
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, 2), atol=1e-5)
    end

    sol = simulate(bg, tspan; u0=[5.0])
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 5, 2), atol=1e-5)
    end

    sol = simulate(bg, tspan; pmap=[1.0, 3.0])
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, 3), atol=1e-5)
    end
end

@testset "Equivalent resistance (DAE)" begin
    R1 = 1.0
    R2 = 2.0
    Req = R1 + R2
    C = 3.0
    τ = Req * C

    r1 = Component(:R, :r1; R=R1)
    r2 = Component(:R, :r2; R=R2)
    c = Component(:C; C=C, q=10)
    kcl = EqualFlow(name=:kcl)
    bg = BondGraph(:RRC)

    add_node!(bg, [c, r1, r2, kcl])
    connect!(bg, c, kcl)
    connect!(bg, kcl, r1)
    connect!(bg, kcl, r2)

    f(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sol = simulate(bg, tspan)
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, τ), atol=1e-5)
    end
end

@testset "π-filter" begin
    Se = Component(:Se, :Pin; es=t->1)

    Pa = EqualEffort(name=:Pa)
    fa = EqualFlow(name=:fa)
    ca = Component(:C, :Ca; C=1, q=1)
    rpa = Component(:R, :Rpa; R=1)

    Pb = EqualEffort(name=:Pb)
    fb = EqualFlow(name=:fb)
    cb = Component(:C, :Cb; C=1, q=2)
    rpb = Component(:R, :Rpb; R=1)

    fs = EqualFlow(name=:fs)
    l = Component(:I, :L; L=1, p=1)
    r = Component(:R, :Rs; R=1)

    rl = Component(:R, :RL; R=1)

    bg = BondGraph(:π_filter)
    add_node!(bg, [Se, Pa, fa, ca, rpa, Pb, fb, cb, rpb, fs, l, r, rl])

    connect!(bg, Se, Pa)
    connect!(bg, Pa, fa)
    connect!(bg, fa, ca)
    connect!(bg, fa, rpa)
    connect!(bg, Pa, fs)
    connect!(bg, fs, l)
    connect!(bg, fs, r)
    connect!(bg, fs, Pb)
    connect!(bg, Pb, fb)
    connect!(bg, fb, cb)
    connect!(bg, fb, rpb)
    connect!(bg, Pb, rl)

    tspan = (0, 100.0)
    sol = simulate(bg, tspan)
    @test sol[1] == [1, 2, 1]
    @test isapprox(sol[end][1], 1.0, atol=1e-5)
    @test isapprox(sol[end][2], 0.5, atol=1e-5)
    @test isapprox(sol[end][3], 0.5, atol=1e-5)
end

@testset "Simulate modular BG" begin
    r = Component(:R; R=1)
    l = Component(:I; L=1, p=1)
    c = Component(:C; C=1, q=1)
    kvl = EqualEffort(name=:kvl)
    SS1 = SourceSensor(name=:SS1)
    SS2 = SourceSensor(name=:SS2)

    bg1 = BondGraph(:RC)
    add_node!(bg1, [r, c, kvl, SS1])
    connect!(bg1, r, kvl)
    connect!(bg1, c, kvl)
    connect!(bg1, SS1, kvl)
    bgn1 = BondGraphNode(bg1)

    bg2 = BondGraph(:L)
    add_node!(bg2, [l, SS2])
    connect!(bg2, l, SS2)
    bgn2 = BondGraphNode(bg2)

    bg = BondGraph()
    add_node!(bg, [bgn1, bgn2])
    connect!(bg, bgn1, bgn2)

    tspan = (0, 10.0)
    sol = simulate(bg, tspan)

    τ = 2
    ω = sqrt(3) / 2
    f(t, τ, ω) = exp(-t / τ) * [
        cos(ω * t) - sqrt(3) * sin(ω * t),
        cos(ω * t) + sqrt(3) * sin(ω * t)
    ]

    for t in [0.0, 0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t), f(t, τ, ω), atol=1e-5)
    end
end

@testset "Driven Filter Circuit" begin
    model = BondGraph("RC")
    C = Component(:C; C=1)
    R = Component(:R; R=1)
    zero_law = EqualEffort()
    C, R, zero_law
    add_node!(model, [C, R, zero_law])
    connect!(model, R, zero_law)
    connect!(model, C, zero_law)

    # Source of flow in the model
    Sf = Component(:Sf)
    add_node!(model, Sf)
    connect!(model, Sf, zero_law)

    # Simulation parameters
    tspan = (0.0, 5.0)
    u0 = [1]

    # Case 1: constant forcing funciton
    Sf.fs = t->3
    ODESystem(model)
    constitutive_relations(model)
    sol = simulate(model, tspan; u0)
    @test isapprox(sol[end], [2.98651], atol=1e-5)

    # Case 2: regular forcing function
    f(t) = sin(2t)
    Sf.fs = f
    sol = simulate(model, tspan; u0)
    @test isapprox(sol[end], [0.23625], atol=1e-5)
end

@testset "Simple Biochemical Simulation" begin
    abc = @reaction_network ABC begin
        1, A + B --> C
    end

    bg_abc = BondGraph(abc)

    sys = ODESystem(bg_abc)
    eqs = constitutive_relations(bg_abc)

    tspan = (0.0, 3.0)
    u0 = [1, 2, 3]
    sol = simulate(bg_abc, tspan; u0)

    @test isapprox(sol[end], [1.23606, 2.23606, 2.76393], atol=1e-5)
end
