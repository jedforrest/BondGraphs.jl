@testset "Set parameters" begin
    re = Component(:Re)

    @test get_default(re, :r) == 1.0
    @test get_default(re, :R) == 8.314
    @test get_default(re, :T) == 310.0

    set_default!(re, :T, 200.0)
    @test get_default(re, :T) == 200.0
end

@testset "Setting non-numeric control variables" begin
    f(t) = sin(2t) # forcing function

    sf = Component(:Sf)
    @test get_default(sf, :fs)(1) ≈ 1

    set_default!(sf, :fs, f)
    @test get_default(sf, :fs) == f
    @test get_default(sf, :fs)(1) ≈ f(1)
end

@testset "Incompatible variables fail" begin
    re = Component(:Re)
    c = Component(:C)
    @test_throws ErrorException set_default!(re, :s, 1.0)
    @test_throws ErrorException set_default!(c, :p, 2.0)
end

@testset "Set initial conditions" begin
    c = Component(:C)
    set_default!(c, :q, 2.0)
    @test get_default(c, :q) == 2.0
end

@testset "Simulate RC circuit" begin
    r = Component(:R)
    c = Component(:C)
    bg = BondGraph(:RC)

    add_node!(bg, [c, r])
    connect!(bg, r, c)

    set_default!(r, :R, 2.0)
    set_default!(c, :C, 1.0)
    set_default!(c, :q, 10.0)

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
    r1 = Component(:R, :r1)
    r2 = Component(:R, :r2)
    c = Component(:C)
    kcl = EqualFlow(name=:kcl)
    bg = BondGraph(:RRC)

    add_node!(bg, [c, r1, r2, kcl])
    connect!(bg, c, kcl)
    connect!(bg, kcl, r1)
    connect!(bg, kcl, r2)

    R1 = 1.0
    R2 = 2.0
    Req = R1 + R2
    C = 3.0
    τ = Req * C

    set_default!(r1, :R, R1)
    set_default!(r2, :R, R2)
    set_default!(c, :C, C)
    set_default!(c, :q, 10.0)

    f(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sol = simulate(bg, tspan)
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, τ), atol=1e-5)
    end
end

@testset "π-filter" begin
    Se = Component(:Se, :Pin)
    set_default!(Se, :es, 1)

    Pa = EqualEffort(name=:Pa)
    fa = EqualFlow(name=:fa)
    ca = Component(:C, :Ca)
    set_default!(ca, :C, 1)
    set_default!(ca, :q, 1)
    rpa = Component(:R, :Rpa)
    set_default!(rpa, :R, 1)

    Pb = EqualEffort(name=:Pb)
    fb = EqualFlow(name=:fb)
    cb = Component(:C, :Cb)
    set_default!(cb, :C, 1)
    set_default!(cb, :q, 2)
    rpb = Component(:R, :Rpb)
    set_default!(rpb, :R, 1)

    fs = EqualFlow(name=:fs)
    l = Component(:I, :L)
    set_default!(l, :L, 1)
    set_default!(l, :p, 1)
    r = Component(:R, :Rs)
    set_default!(r, :R, 1)

    rl = Component(:R, :RL)
    set_default!(rl, :R, 1)

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
    r = Component(:R)
    set_default!(r, :R, 1)
    l = Component(:I)
    set_default!(l, :L, 1)
    set_default!(l, :p, 1)
    c = Component(:C)
    set_default!(c, :C, 1)
    set_default!(c, :q, 1)
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

    constitutive_relations(bg)

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
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()
    C, R, zero_law
    add_node!(model, [C, R, zero_law])
    connect!(model, R, zero_law)
    connect!(model, C, zero_law)
    set_default!(C, :C, 1.0)
    set_default!(R, :R, 1.0)

    # Source of flow in the model
    Sf = Component(:Sf)
    add_node!(model, Sf)
    connect!(model, Sf, zero_law)

    # Simulation parameters
    tspan = (0.0, 5.0)
    u0 = [1]

    # Case 1: constant forcing funciton
    set_default!(Sf, :fs, -3)
    ODESystem(model)
    constitutive_relations(model)
    sol = simulate(model, tspan; u0)
    @test isapprox(sol[end], [2.98651], atol=1e-5)

    # Case 2: regular forcing function
    f(t) = -sin(2t)
    set_default!(Sf, :fs, f)
    sol = simulate(model, tspan; u0)
    @test isapprox(sol[end], [0.23625], atol=1e-5)

    # @register_symbolic needs to be called at the top level
    # Therefore this test fails, despite the code working in REPL

    # Case 3: arbitrary forcing function
    # g(t) = t % 1 <= 0.5 ? 2 : 0
    # @register_symbolic g(t)
    # set_default!(Sf, :fs, g)
    # sol = simulate(model, tspan; u0)
    # @test isapprox(sol[end], [-0.70709], atol=1e-5)
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