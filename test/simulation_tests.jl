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
