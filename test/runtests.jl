using BondGraphs
using LightGraphs
using Test

@testset "BondGraphs.jl" begin
    c1 = Component(:C)
    @test c1.name == "C"

    c2 = Component(:C, "newC", 1)
    @test c1 != c2

    j = Junction(:J)
    b = Bond(c1, j)
    @test src(b) == c1
    @test dst(b) == j

    bg = BondGraph()
    @test bg.metamodel == :BG
    @test isempty(bg.nodes)
end
