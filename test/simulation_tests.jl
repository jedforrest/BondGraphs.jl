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