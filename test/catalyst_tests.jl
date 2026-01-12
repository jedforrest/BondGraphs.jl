@testset "Chemical reaction" begin
    model = BondGraph(:Chemical)
    A = Component(:C, :A)
    B = Component(:C, :B)
    C = Component(:C, :C)
    D = Component(:C, :D)
    Re = Component(:Re, :Reaction, numports = 2)
    J_AB = EqualFlow()
    J_CD = EqualFlow()

    add_node!(model, [A, B, C, D, Re, J_AB, J_CD])
    connect!(model, A, J_AB)
    connect!(model, B, J_AB)
    connect!(model, C, J_CD)
    connect!(model, D, J_CD)

    # Connecting junctions to specific ports in Re
    connect!(model, (Re, 2), J_CD)
    @test ports(Re) == Dict(1 => false, 2 => true)
    connect!(model, J_AB, (Re, 1))
    @test ports(Re) == Dict(1 => true, 2 => true)

    @test nv(model) == 7
    @test ne(model) == 6
end

# FIXME just for chemical species really
@testset "Merging components" begin
    bg = RCI()
    C = bg.C
    R = bg.R

    newC = Component(:C, :newC)
    newR = Component(:R, :newR)
    add_node!(bg, [newC, newR])
    connect!(bg, newC, newR)

    merge_nodes!(bg, C, newC)
    @test isempty(getnodes(bg, "C:newC"))

    merge_nodes!(bg, R, newR; junction = EqualFlow())
    @test length(getnodes(bg, EqualFlow)) == 1
    @test length(getnodes(bg, EqualEffort)) == 2
    @test nv(bg) == 7
    @test ne(bg) == 7
end

@testset "Simple Reaction System" begin
    rn = @reaction_network ABC begin
        1, A + B --> C
    end

    bg_rn = BondGraph(rn)

    @test bg_rn.name == "ABC"
    @test nv(bg_rn) == 5
    @test ne(bg_rn) == 4

    # Checking that a "reverse" direction bond is included
    bond_rev = Bond((bg_rn.C, 1), (bg_rn.R1, 2))
    @test bond_rev in bonds(bg_rn)

    @test Graphs.degree(bg_rn) == [2, 3, 1, 1, 1]
end

@testset "Reversible MM" begin
    rn = @reaction_network MM_reversible begin
        (1, 1), E + S <--> C
        (1, 1), C <--> E + P
    end

    bg_rn = BondGraph(rn; chemostats = ["S", "P"])

    @test name(bg_rn) == "MM_reversible"
    @test nv(bg_rn) == 10
    @test ne(bg_rn) == 10

    @test Graphs.degree(bg_rn) == [2, 3, 1, 1, 1, 2, 3, 1, 3, 3]
end

@testset "Stoichiometry Test" begin
    rn = @reaction_network Stoichiometry begin
        1, 3A + 2B --> 5C
    end

    bg_rn = BondGraph(rn)

    @test nv(bg_rn) == 8
    @test ne(bg_rn) == 7

    tfs = filter(n -> type(n) == "TF", bg_rn.nodes)
    @test repr.(tfs) == ["TF:tf1", "TF:tf2", "TF:tf3"]
end

@testset "SERCA" begin
    rn = @reaction_network SERCA begin
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
    bg_rn = BondGraph(rn; chemostats)

    @test nv(bg_rn) == 46
    @test ne(bg_rn) == 49
end
