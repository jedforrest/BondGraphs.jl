@testset "Set parameters" begin
    re = Component(:Re)

    @test get_parameter(re, :r) == 1.0
    @test get_parameter(re, :R) == 8.314
    @test get_parameter(re, :T) == 310.0

    set_parameter!(re, :T, 200.0)
    @test get_parameter(re, :T) == 200.0
end

@testset "Incompatible parameter fail" begin
    re = Component(:Re)
    @test_throws ErrorException set_parameter!(re, :s, 1.0)
end

@testset "Set initial conditions" begin
    c = Component(:C)
    set_initial_value!(c, :q, 2.0)
    @test get_initial_value(c, :q) == 2.0
end

@testset "Missing state variable fail" begin
    c = Component(:C)
    @test_throws ErrorException set_parameter!(c, :q, 2.0)
    @test_throws ErrorException set_initial_value!(c, :p, 2.0)
end

@testset "Simulate RC circuit" begin
    r = Component(:R)
    c = Component(:C)
    bg = BondGraph(:RC)

    add_node!(bg, [c, r])
    connect!(bg, r, c)

    set_parameter!(r, :R, 2.0)
    set_parameter!(c, :C, 1.0)
    set_initial_value!(c, :q, 10.0)

    f(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sol = simulate(bg, tspan)
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, 2), atol = 1e-5)
    end

    sol = simulate(bg, tspan; u0 = [5.0])
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 5, 2), atol = 1e-5)
    end

    sol = simulate(bg, tspan; pmap = [1.0, 3.0])
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, 3), atol = 1e-5)
    end
end

@testset "Equivalent resistance (DAE)" begin
    r1 = Component(:R, :r1)
    r2 = Component(:R, :r2)
    c = Component(:C)
    kcl = EqualFlow(name = :kcl)
    bg = BondGraph(:RRC)
    
    add_node!(bg, [c, r1, r2, kcl])
    connect!(bg, c, kcl)
    connect!(bg, kcl, r1)
    connect!(bg, kcl, r2)
    
    R1 = 1.0
    R2 = 2.0
    Req = R1+R2
    C = 3.0
    τ = Req*C
    
    set_parameter!(r1, :R, R1)
    set_parameter!(r2, :R, R2)
    set_parameter!(c, :C, C)
    set_initial_value!(c, :q, 10.0)

    f(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sol = simulate(bg, tspan)
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f(t, 10, Req*C), atol = 1e-5)
    end
end

@testset "π-filter" begin
    Se = Component(:Se, :Pin); set_parameter!(Se, :e, 1)

    Pa = EqualEffort(name = :Pa)
    fa = EqualFlow(name = :fa)
    ca = Component(:C, :Ca); set_parameter!(ca, :C, 1); set_initial_value!(ca, :q, 1)
    rpa = Component(:R, :Rpa); set_parameter!(rpa, :R, 1)
    
    Pb = EqualEffort(name = :Pb)
    fb = EqualFlow(name = :fb)
    cb = Component(:C, :Cb); set_parameter!(cb, :C, 1); set_initial_value!(cb, :q, 2)
    rpb = Component(:R, :Rpb); set_parameter!(rpb, :R, 1)
    
    fs = EqualFlow(name = :fs)
    l = Component(:I, :L); set_parameter!(l, :L, 1); set_initial_value!(l, :p, 1)
    r = Component(:R, :Rs); set_parameter!(r, :R, 1)
    
    rl = Component(:R, :RL); set_parameter!(rl, :R, 1)
    
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
    
    tspan = (0,100.0)
    sol = simulate(bg,tspan)
    @test sol[1] == [1,2,1]
    @test isapprox(sol[end][1], 1.0, atol=1e-5)
    @test isapprox(sol[end][2], 0.5, atol=1e-5)
    @test isapprox(sol[end][3], 0.5, atol=1e-5)
end

@testset "Simulate modular BG" begin
    r = Component(:R); set_parameter!(r, :R, 1)
    l = Component(:I); set_parameter!(l, :L, 1); set_initial_value!(l, :p, 1)
    c = Component(:C); set_parameter!(c, :C, 1); set_initial_value!(c, :q, 1)
    kvl = EqualEffort(name = :kvl)
    SS1 = SourceSensor(name = :SS1)
    SS2 = SourceSensor(name = :SS2)
    
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
    
    tspan = (0,10.0)
    sol = simulate(bg,tspan)
    
    τ = 2
    ω = sqrt(3)/2
    f(t,τ,ω) = exp(-t/τ)*[
        cos(ω*t) - sqrt(3)*sin(ω*t),
        cos(ω*t) + sqrt(3)*sin(ω*t)
    ]
    
    for t in [0.0, 0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t), f(t,τ,ω), atol = 1e-5)
    end
end