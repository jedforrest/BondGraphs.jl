@testset "Chemical reaction" begin
    @named X = chemicalspecies()
    @named Y = chemicalspecies()
    @named Z = chemicalspecies()
    @named W = chemicalspecies()
    @named Re = reaction()
    J_XY = EqualFlow(name=:j1)
    J_ZW = EqualFlow(name=:j2)

    model = BondGraph([X, Y, Z, W, Re, J_XY, J_ZW], name=:Chemical)
    connect!(model, X, J_XY)
    connect!(model, Y, J_XY)
    connect!(model, Z, J_ZW)
    connect!(model, W, J_ZW)

    # Connecting junctions to specific ports in Re
    connect!(model, Re[2], J_ZW)
    @test is_connected.(ports(Re)) == [false, true]
    connect!(model, J_XY, Re[1])
    @test is_connected.(ports(Re)) == [true, true]

    @test nv(model) == 7
    @test ne(model) == 6
end

@testset "Merging components" begin
    bg = RCI()
    C = bg[:c]
    R = bg[:r]

    @named newC = capacitor()
    @named newR = resistor()

    add_comp!(bg, newC)
    add_comp!(bg, newR)
    connect!(bg, newC, newR)

    merge_comps!(bg, C, newC)
    @test !(newC in bg)

    merge_comps!(bg, R, newR; junctiontype = EqualFlow)
    @test length(junctions(bg)) == 3
    @test nv(bg) == 7
    @test ne(bg) == 7
end

# TODO CONTINUE FROM HERE
# rewrite catalyst -> bondgraph function
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
