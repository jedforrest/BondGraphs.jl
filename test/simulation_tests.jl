@testset "Set parameters" begin
    set_library!(biochemical_library)
    
    re = Component(:Re)
    @parameters r, R, T
    set_param!(re,r,1.0)
    set_param!(re,R,8.314)
    set_param!(re,T,310.0)

    @test default_value(re,r) == 1.0
    @test default_value(re,R) == 8.314
    @test default_value(re,T) == 310.0

    set_library!()
end

@testset "Incompatible parameter fail" begin
    re = Component(:Re)
    @parameters s
    @test_throws ErrorException set_param!(re,s,1.0)
end

@testset "Set initial conditions" begin
    c = Component(:C)
    @variables q(t)
    set_initial_value!(c,q,2.0)
    @test default_value(c,q) == 2.0
end

@testset "Missing state variable fail" begin
    c = Component(:C)
    @parameters q
    @variables p(t)
    @test_throws ErrorException set_initial_value!(c,q,2.0)
    @test_throws ErrorException set_initial_value!(c,p,2.0)
end

@testset "Simulate RC circuit" begin
    r = Component(:R)
    c = Component(:C)
    bg = BondGraph(:RC)

    add_node!(bg, [c, r])
    connect!(bg, r, c)

    @parameters R, C
    @variables q(t)
    set_param!(r,R,2.0)
    set_param!(c,C,1.0)
    set_initial_value!(c,q,10.0)

    f(x,a,τ) = a*exp(-x/τ)

    tspan = (0.0,10.0)
    sol = simulate(bg,tspan)
    for t in [0.5,1.0,5.0,10.0]
        @test isapprox(sol(t)[1], f(t,10,2), atol=1e-5)
    end

    sol = simulate(bg,tspan; u0=[5.0])
    for t in [0.5,1.0,5.0,10.0]
        @test isapprox(sol(t)[1], f(t,5,2), atol=1e-5)
    end

    sol = simulate(bg,tspan; pmap=[1.0,3.0])
    for t in [0.5,1.0,5.0,10.0]
        @test isapprox(sol(t)[1], f(t,10,3), atol=1e-5)
    end
end