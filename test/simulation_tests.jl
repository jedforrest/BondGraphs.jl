@testset "Set parameters" begin
    re = new(:Re)
    @parameters r, R, T
    set_param!(re,r,1.0)
    set_param!(re,R,8.314)
    set_param!(re,T,310.0)

    @test default_value(re,r) == 1.0
    @test default_value(re,R) == 8.314
    @test default_value(re,T) == 310.0
end

@testset "Incompatible parameter fail" begin
    re = new(:Re)
    @parameters s
    @test_throws ErrorException set_param!(re,s,1.0)
end

@testset "Set initial conditions" begin
    c = new(:C)
    @variables q(t)
    set_initial_value!(c,q,2.0)
    @test default_value(c,q) == 2.0
end

@testset "Missing state variable fail" begin
    c = new(:C)
    @parameters q
    @variables p(t)
    @test_throws ErrorException set_initial_value!(c,q,2.0)
    @test_throws ErrorException set_initial_value!(c,p,2.0)
end

@testset "Simulate RC circuit" begin
    r = new(:R)
    c = new(:C)
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

# Notes: The SciMLBase package contains a getindex function for syms.
# Within this function, sym_to_index(sym,A) will return the index of the variable.
@testset "Extract specific variables" begin
    # Make a model of the reaction A + B ⇌ C + D
    C_A = new(:ce,:A)
    C_B = new(:ce,:B)
    C_C = new(:ce,:C)
    C_D = new(:ce,:D)
    re = new(:re,:r)
    AB = EqualFlow(name=:AB)
    CD = EqualFlow(name=:CD)
    
    bg = BondGraph()
    add_node!(bg,[C_A,C_B,C_C,C_D,re,AB,CD])
    connect!(bg,C_A,AB)
    connect!(bg,C_B,AB)
    connect!(bg,AB,re; dstportindex=1)
    connect!(bg,re,CD; srcportindex=2)
    connect!(bg,CD,C_C)
    connect!(bg,CD,C_D)
    
    @parameters t, r, k
    @variables q(t)
    set_param!(C_A,k,1.0)
    set_param!(C_B,k,2.0)
    set_param!(C_C,k,1.0)
    set_param!(C_D,k,2.0)
    set_param!(re,r,1.0)
    set_initial_value!(C_A,q,1.0)
    set_initial_value!(C_B,q,4.0)
    set_initial_value!(C_C,q,2.0)
    set_initial_value!(C_D,q,1.0)
    
    tspan = (t0,t1) = (0.0,10.0)
    sol = simulate(bg,tspan)

    xA = sol["C:A"]
    @test xA[1] = 1.0
    @test xA[end] ≈ 0.75
    
    # Test that system is near equilibrium
    @test sol("C:A";t=t1)*sol("C:B";t=t1)/sol("C:C";t=t1)/sol("C:D";t=t1) ≈ 1
end