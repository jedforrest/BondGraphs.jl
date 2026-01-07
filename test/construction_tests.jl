using Test
using BondGraphs
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables e(t) f(t) p(t) q(t)
@parameters R C L

@testset "power variables" begin
    using ModelingToolkit: get_connection_type

    e, f, _ = power_variables()
    @test get_connection_type(e) == Effort
    @test get_connection_type(f) == Flow

    vars = power_variables(e="F", f="v", p="p", q="x")
    @test tosymbol.(vars, escape=false) == [:F, :v, :p, :x]
end

@testset "AbstractElements" #= setup=[Setup] =# begin

    r_element = DissipatorElement([e ~ R * f], [e], [f])
    c_element = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q])
    i_element = DynamicStorageElement([D(p) ~ e, p ~ L * f], [e], [f], [p])

    @test isequal(r_element.efforts, [e])
    @test isequal(c_element.states, [q])
    @test isequal(equations(i_element.sys), [D(p) ~ e, p ~ L * f])
end

@testset "Creating Components" #= setup=[Setup] =# begin
    @variables e(t) f(t) p(t) q(t)
    @parameters R C L

    c_element = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q])
    r_element = DissipatorElement([e ~ R * f], [e], [f])

    @named c_comp = Component(c_element; q = 2)
    @named r_comp = Component(r_element; R = 5)

    @test elementtype(c_comp) == StaticStorageElement
    @test name(c_comp) == :c_comp
    @test repr(c_comp) == "C::c_comp"

    @test isequal(collect(defaults(c_comp.sys)), [q => 2])
    @test isequal(collect(defaults(r_comp.sys)), [R => 5])
end

# TODO CONTINUE FROM HERE
@testset "Creating Junctions" #= setup=[Setup] =# begin
    EqE_1 = EqualEffort()
    EqE_2 = EqualEffort(name = "foo")
    EqF = EqualFlow()

    @test name(EqE_1) == "𝟎"
    @test name(EqE_2) == "foo"
    @test name(EqF) == "𝟏"
end

@testset "BondGraph Construction" #= setup=[Setup] =# begin
    model = BondGraph(:RC)
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    @test R in model.nodes
    @test C in model.nodes
    @test zero_law in model.nodes

    b1 = connect!(model, R, zero_law)
    b2 = connect!(model, zero_law, C)
    @test b1 in model.bonds
    @test b2 in model.bonds
end

@testset "Graph construction" #= setup=[Setup] =# begin
    c1 = Component(:C)
    c2 = Component(:R)
    c3 = Component(:I)
    j = EqualFlow()
    bg = BondGraph()

    add_vertex!(bg, c1)
    add_vertex!(bg, c2)
    add_vertex!(bg, c3)
    add_vertex!(bg, j)
    add_edge!(bg, (c1, 1), (j, 1))
    add_edge!(bg, (j, 1), (c2, 1)) # junction index is '1' here as a quick fix
    add_edge!(bg, (j, 1), (c3, 1)) # junction index is '1' here as a quick fix

    # adding components and bonds
    @test ne(bg) == 3
    @test nv(bg) == 4
    @test has_vertex(bg, c1)
    @test has_edge(bg, vertex(j), vertex(c2))

    # example graph functions
    @test Δ(bg) == 3
    @test Graphs.density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]

    # removing components and bonds
    rem_edge!(bg, c3, j)
    rem_vertex!(bg, c3)
    @test ne(bg) == 2
    @test nv(bg) == 3
end

@testset "BondGraph Modification" #= setup=[Setup] =# begin
    model = BondGraph(:RCI)
    C = Component(:C)
    R = Component(:R)
    I = Component(:I)
    SS = Component(:SS)
    zero_law = EqualEffort()
    one_law = EqualFlow()

    add_node!(model, [C, R, I, SS, zero_law, one_law])
    remove_node!(model, [SS, one_law])
    @test !(SS in model.nodes)
    @test !(one_law in model.nodes)

    connect!(model, R, zero_law)
    connect!(model, C, zero_law)

    I_port_info = BondGraphs.port_info(I)
    @test I_port_info == (I, 1)
    @test I_port_info == BondGraphs.port_info(I_port_info)

    @test ports(I) == Dict(1 => false)
    b1 = connect!(model, zero_law, I)
    @test ne(model) == 3
    @test b1 in bonds(model)
    @test ports(I) == Dict(1 => true)

    disconnect!(model, I, zero_law) # tests disconnect when node order is swapped
    @test ne(model) == 2
    @test !(b1 in bonds(model))
    @test I.ports == Dict(1 => false)

    connect!(model, zero_law, I)
    swap!(model, zero_law, one_law)
    @test ports(I) == Dict(1 => true)
    @test one_law in nodes(model)
    @test inneighbors(model, one_law) == [R, C]
    @test outneighbors(model, one_law) == [I]
end

@testset "Construction Failure" #= setup=[Setup] =# begin
    model = BondGraph(:RC)
    C = Component(:C)
    R = Component(:R)
    zero_law = EqualEffort()

    add_node!(model, [R, C, zero_law])
    @test_logs (:warn, "Node 'R' already in model") add_node!(model, R)
    @test_logs (:warn, "Node '𝟎_3' already in model") add_node!(model, zero_law)

    bond = connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, R, zero_law)
    @test_throws ErrorException connect!(model, C, R)

    one_law = EqualFlow()
    @test_logs (:warn, "Node '𝟏' not in model") remove_node!(model, one_law)

    tf = Component(:TF)
    add_node!(model, tf)
    @test_throws ErrorException swap!(model, tf, C)

    # if inserting a node fails, the original nodes should still remain connected
    @test has_edge(model, bond)
    @test_throws ErrorException insert_node!(model, bond, Component(:I))
    @test has_edge(model, bond)
end

@testset "Chemical reaction" #= setup=[Setup] =# begin
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

@testset "Inserting Nodes" #= setup=[Setup] =# begin
    bg = RCI()

    c, r, J0 = bg.nodes[[1, 2, 5]]

    bondc0 = getbonds(bg, c, J0)[1]
    bondr0 = getbonds(bg, r, J0)[1]

    tf = Component(:TF, numports = 2)
    insert_node!(bg, bondc0, tf)
    insert_node!(bg, bondr0, EqualFlow())

    @test tf in bg.nodes
    @test nv(bg) == 7
    @test ne(bg) == 6
end

@testset "Merging components" #= setup=[Setup] =# begin
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

@testset "Simplifying Junctions" #= setup=[Setup] =# begin
    bg = RCI()
    C, R, I, SS, J0 = bg.nodes

    J0_new_1 = EqualEffort(; name = :new0_1)
    J0_new_2 = EqualEffort(; name = :new0_2)
    insert_node!(bg, (C, J0), J0_new_1)
    insert_node!(bg, (R, J0), J0_new_2)
    connect!(bg, J0_new_1, J0_new_2)

    J1_new_1 = EqualFlow(; name = :new1_1)
    J1_new_2 = EqualFlow(; name = :new1_2)
    add_node!(bg, J1_new_1)
    connect!(bg, J0, J1_new_1)
    insert_node!(bg, (SS, J0), J1_new_2)

    # Removing junction redundancies
    @test length(getnodes(bg, EqualFlow)) == 2
    simplify_junctions!(bg, squash_identical = false)
    @test length(getnodes(bg, EqualFlow)) == 0
    @test nv(bg) == 7
    @test ne(bg) == 7

    # Squashing junction duplicates into a single junction
    simplify_junctions!(bg)
    @test length(getnodes(bg, EqualEffort)) == 1
    @test nv(bg) == 5
    @test ne(bg) == 4
end

@testset "BondGraphNodes" #= setup=[Setup] =# begin
    C = Component(:C, "C")
    bg1 = BondGraph("first")
    bg2 = BondGraph("second")
    bg3 = BondGraph("third")
    main = BondGraph("Main")

    bgn1 = BondGraphNode(bg1)
    bgn2 = BondGraphNode(bg2)
    bgn3 = BondGraphNode(bg3)

    @test bgn1.bondgraph === bg1
    @test bgn1.type === "BG"
    @test bgn1.name === bg1.name
    @test bgn1.ports == Dict()

    add_node!(bg1, C)
    add_node!(bg2, bgn1)
    add_node!(bg3, bgn2)
    add_node!(main, bgn3)

    @test main.third.second.first.C === C

    C2 = Component(:C, "C") # Same name
    add_node!(bg1, C2)
    @test main.third.second.first.C == [C, C2]
end

@testset "Conversion to Other Graphs" #= setup=[Setup] =# begin
    bg = RCI()
    g = SimpleGraph(bg)
    dg = SimpleDiGraph(bg)

    bg_adj, g_adj, dg_adj = adjacency_matrix.([bg, g, dg])
    @test dg_adj == bg_adj
    @test g_adj == bg_adj + bg_adj' # A + A' forms undirected graph adj matrix
end
