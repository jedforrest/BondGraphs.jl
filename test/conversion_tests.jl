@testset "Reaction Systems" begin
    rn = @reaction_network ABC begin
        1, A + B --> C
    end

    bg_rn = BondGraph(rn)

    @test bg_rn.name == "ABC"
    @test nv(bg_rn) == 5
    @test ne(bg_rn) == 4
    
    @test any(n -> n.name == "R1", bg_rn.nodes)
    @test any(n -> typeof(n) == Junction && n.metamodel == :𝟏, bg_rn.nodes)

    @test length(getnodes(bg_rn, :Ce)) == 3
    @test length(getnodes(bg_rn, :𝟏)) == 1
    @test length(getnodes(bg_rn, :Re)) == 1

    # LightGraphs degree function
    @test degree(bg_rn) == [2, 3, 1, 1, 1]
end