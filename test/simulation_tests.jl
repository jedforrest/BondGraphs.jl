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
    Se = Component(:Se, :Pin; es=t -> 1)

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
    Sf.fs = t -> 3
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

    tspan = (0.0, 3.0)
    u0 = [1, 2, 3]
    sol = simulate(bg_abc, tspan; u0)

    @test isapprox(sol[end], [1.23606, 2.23606, 2.76393], atol=1e-5)
end

@testset "Stoichiometry Simulation" begin
    rn = @reaction_network A2B begin
        1, A --> 2B
    end
    bg = BondGraph(rn)
    sol = simulate(bg, (0.0, 1.0); u0=[1, 0])
    @test isapprox(sol[end], [0.61969, 0.76062], atol=1e-5)
end

@testset "Reversible Michaelis-Menten" begin
    rn_mm = @reaction_network MM_reversible begin
        (1, 1), E + S <--> C
        (1, 1), C <--> E + P
    end
    bg_mm = BondGraph(rn_mm; chemostats=["S", "P"])

    bg_mm.S.xs = t -> 2

    sol = simulate(bg_mm, (0.0, 3.0); u0=[1, 2])

    @test isapprox(sol[end], [1.2, 1.8], atol=1e-5)
end

@testset "SERCA (stiff equations)" begin
    rn_serca = @reaction_network SERCA begin
        (1, 1), P1 + MgATP <--> P2
        (1, 1), P2 + H <--> P2a
        (1, 1), P2 + 2Cai <--> P4
        (1, 1), P4 <--> P5 + 2H
        (1, 1), P5 <--> P6 + MgADP
        (1, 1), P6 <--> P8 + 2Casr
        (1, 1), P8 + 2H <--> P9
        (1, 1), P9 <--> P10 + H
        (1, 1), P10 <--> P1 + Pi
    end

    chemostats = ["MgATP", "MgADP", "Pi", "H", "Cai", "Casr"]
    bg_serca = BondGraph(rn_serca; chemostats)


    reaction_rates = [
        :R1 => 0.00053004,
        :R2 => 8326784.0537,
        :R3 => 1567.7476,
        :R4 => 1567.7476,
        :R5 => 3063.4006,
        :R6 => 130852.3839,
        :R7 => 11612934.8748,
        :R8 => 11612934.8748,
        :R9 => 0.049926
    ]
    for (reaction, rate) in reaction_rates
        getproperty(bg_serca, reaction).r = rate
    end

    species_affinities = [
        :P1 => 5263.6085,
        :P2 => 3803.6518,
        :P2a => 3110.4445,
        :P4 => 16520516.1239,
        :P5 => 0.82914,
        :P6 => 993148.433,
        :P8 => 37.7379,
        :P9 => 2230.2717,
        :P10 => 410.6048,
        :Cai => 1.9058,
        :Casr => 31.764,
        :MgATP => 244.3021,
        :MgADP => 5.8126e-7,
        :Pi => 0.014921,
        :H => 1862.5406
    ]
    for (species, affinity) in species_affinities
        getproperty(bg_serca, species).K = affinity
    end

    chemostat_amounts = [
        :Cai => t -> 0.0057,
        :Casr => t -> (0.05 + 0.01t)*2.28,
        :H => t -> 0.004028,
        :MgADP => t -> 1.3794,
        :MgATP => t -> 3.8,
        :Pi => t -> 570
    ]
    for (chemostat, amount) in chemostat_amounts
        getproperty(bg_serca, chemostat).xs = amount
    end

    initial_conditions = [
        :P1 => 0.000483061870385487,
        :P2 => 0.0574915174273067,
        :P2a => 0.527445119834607,
        :P4 => 1.51818391164022e-09,
        :P5 => 0.000521923287622898,
        :P6 => 7.80721128535043e-05,
        :P8 => 0.156693953834181,
        :P9 => 0.149232225342376,
        :P10 => 0.108044124948978
    ]
    for (species, ic) in initial_conditions
        getproperty(bg_serca, species).q = ic
    end

    tspan = (0., 200.)
    sol = simulate(bg_serca, tspan; solver=Rosenbrock23())

    # calculated using the same code, verified by plot from BGT tutorial
    real_solution = [
        4.4404656222265794e-5,
        0.09777422826977565,
        0.8970112784324162,
        2.6596475539704174e-9,
        0.0009426424413096248,
        0.001015195974904865,
        0.001212098675876874,
        0.0011543788312157496,
        0.0008357702283899367,
    ]
    @test isapprox(sol[end], real_solution, atol=1e-5)
end
