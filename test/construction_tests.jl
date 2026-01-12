using Test
using BondGraphs
using BondGraphs: is_connected
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Graphs, MetaGraphsNext

@testset "power variables" begin
    e, f, _ = power_variables()
    @test ModelingToolkit.get_connection_type(e) == Effort
    @test ModelingToolkit.get_connection_type(f) == Flow

    vars = power_variables(e="F", f="v", p="p", q="x")
    @test tosymbol.(vars, escape=false) == [:F, :v, :p, :x]
end

@testset "AbstractElements" begin
    @variables e(t) f(t) p(t) q(t)
    @parameters R C L

    r_element = DissipatorElement([e ~ R * f], [e], [f])
    c_element = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q])
    i_element = DynamicStorageElement([D(p) ~ e, p ~ L * f], [e], [f], [p])

    @test isequal(efforts(r_element), [e])
    @test isequal(flows(c_element), [f])
    @test isequal(states(c_element), [q])
    @test isequal(constitutive_relations(i_element), [D(p) ~ e, p ~ L * f])
end

@testset "Ports" begin
    port = Port(:port_1, :parent, e=:V, f=:I)

    @test nameof(system(port)) == :parent₊port_1
    @test nameof(system(port, namespaced=false)) == :port_1

    @test is_connected(port) == false
    connect!(port)
    @test is_connected(port) == true

    @test repr(effort(port)) == "port_1₊V(t)"
    @test repr(flow(port)) == "port_1₊I(t)"
end

@testset "Creating Components" begin
    @variables e(t) f(t) q(t)
    @parameters C
    eqs = [D(q) ~ f, q ~ C * e]

    @named capacitor = StaticStorageElement(eqs, [e], [f], [q]; q = 2, C = 5)

    @test typeof(capacitor) == StaticStorageElement
    @test name(capacitor) == :capacitor
    @test repr(capacitor) == "C::capacitor"

    @test capacitor.q == 2
    @test capacitor.C == 5

    port = capacitor.ports[1]
    @test capacitor[1] == port
    @test capacitor[:port_1] == port
end

@testset "Creating Junctions" begin
    j0 = EqualEffort()
    j1 = EqualFlow()
    @test repr(j0) == "EqualEffort"
    @test repr(j1) == "EqualFlow"

    @test name(j0) == :𝟎
    @test name(j1) == :𝟏

    @test length(ports(j0)) == length(ports(j1)) == 1

    push!(j0.ports, Port(:port_2, :j0))
    push!(j0.ports, Port(:port_3, :j0))

    cr = constitutive_relations(j0)
    @test repr(cr[1]) == "port_1₊e(t) ~ port_2₊e(t)"
    @test repr(cr[2]) == "port_1₊e(t) ~ port_3₊e(t)"
    @test repr(cr[3]) == "port_3₊f(t) + port_1₊f(t) + port_2₊f(t) ~ 0"
end

@testset "BondGraph Construction" begin
    # imported from library
    using BondGraphs: resistor, capacitor, KCL
    @named rcomp = resistor()
    @named ccomp = capacitor()
    @named kcl = KCL()
    b1 = Bond(rcomp, kcl)
    b2 = Bond(ccomp, kcl)
    bg = BondGraph([rcomp, ccomp, kcl], [b1, b2], name="RC Circuit")

    @test name(bg) == Symbol("RC Circuit")
    @test components(bg) == [rcomp, ccomp, kcl]
    @test elements(bg) == [rcomp, ccomp]
    @test junctions(bg) == [kcl]
    @test bonds(bg) == [b1, b2]

    # compare to bonds only construction
    # bg2 = BondGraph([b1, b2], name="RC Circuit")
    # TODO
end

@testset "Graph functions" begin
    using BondGraphs: resistor, capacitor, inductor, KVL
    @named r = resistor()
    @named c = capacitor()
    @named i = inductor()
    @named kvl = KVL()
    b1 = Bond(c, kvl)
    b2 = Bond(kvl, r)
    b3 = Bond(kvl, i)
    bg = BondGraph([c, r, i, kvl], [b1, b2, b3], name=:RCI)

    @test eltype(BondGraph) == Int
    @test eltype(bg) == Int
    @test edgetype(BondGraph) == Graphs.SimpleEdge{Int}
    @test edgetype(bg) == Graphs.SimpleEdge{Int}

    @test is_directed(bg)

    @test collect(labels(bg.graph)) == [:c, :r, :i, :kvl]
    @test collect(edge_labels(bg.graph)) == [(:c, :kvl), (:kvl, :r), (:kvl, :i)]

    @test ne(bg) == 3
    @test nv(bg) == 4

    # example graph functions
    @test Δ(bg) == 3
    @test Graphs.density(bg) == 0.25
    @test Array(adjacency_matrix(bg)) == [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 1 1 0]

    # removing components and bonds
    rem_edge!(bg.graph, 1, 4)  # c, kvl
    rem_vertex!(bg.graph, 1)  # c
    @test ne(bg) == 2
    @test nv(bg) == 3

    # metagraph getindex
    @test bg[] == "RCI"
    @test bg[:r] == r
    @test bg[:kvl, :i] == b3
end

# TODO CONTINUE FROM HERE
@testset "BondGraph Modification" begin
    using BondGraphs: resistor, capacitor, inductor, voltagesource, KCL, KVL
    @named r = resistor()
    @named c = capacitor()
    @named i = inductor()
    @named v = voltagesource()
    @named kcl = KCL()  # 0-junction
    @named kvl = KVL()  # 1-junction
    @named model = BondGraph()

    model[:r] = r
    model[:c] = c
    model[:i] = i
    for newcomp in [v, kcl, kvl]
        add_comp!(model, newcomp)
    end
    @test nv(model) == 6

    model[:kvl, :kcl] = Bond(kvl, kcl)
    @test length(kcl.ports) == 1 && is_connected(kcl[1])

    @test remove_comp!(model, kvl)
    @test length(kcl.ports) == 1 && !is_connected(kcl[1])
    @test nv(model) == 5
    @test ne(model) == 0

    connect!(model, r, kcl)
    connect!(model, c, kcl)
    connect!(model, kcl, i)
    @test ne(model) == 3
    @test length(kcl.ports) == 3 && all(is_connected, kcl.ports)

    disconnect!(model, kcl, i)
    @test ne(model) == 2
    @test length(kcl.ports) == 3 && !is_connected(kcl[3])
    connect!(model, kcl, i)

    swap!(model, kcl, kvl)
    @test (kvl in components(model)) && !(kcl in components(model))
    @test ne(model) == 3
    @test inneighbor_comps(model, kvl) == [:r, :c]
    @test outneighbor_comps(model, kvl) == [:i]
end

# TODO move to runtests.jl
using BondGraphs: resistor, capacitor, inductor, voltagesource, KVL
function RCI()
    @named r = resistor()
    @named c = capacitor()
    @named i = inductor()
    @named v = voltagesource()
    @named kvl = KVL()
    b2 = Bond(r, kvl)
    b1 = Bond(c, kvl)
    b3 = Bond(kvl, i)
    b4 = Bond(kvl, v)
    BondGraph([c, r, i, v, kvl], [b1, b2, b3, b4], name=:RCI)
end

### TODO CONTINUE FROM HERE
@testset "Inserting Nodes" begin
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

@testset "Construction Failure" begin
    using BondGraphs: resistor, capacitor, KVL
    @named rcomp = resistor()
    @named ccomp = capacitor()
    @named kvl = KVL()
    b1 = Bond(rcomp, kvl)
    b2 = Bond(ccomp, kvl)
    bg = BondGraph([rcomp, ccomp, kvl], [b1, b2], name="RC")

    @test_logs (:warn, "Component 'rcomp' already in model") add_comp!(bg, rcomp)

    @named v = voltagesource()
    @test_throws ErrorException connect!(bg, v, kvl)
    @test_throws ErrorException connect!(bg, ccomp, rcomp)

    @test_logs (:warn, "Component 'v' not in model") remove_comp!(bg, v)

    # FIXME
    # tf = Component(:TF)
    # add_node!(model, tf)
    # @test_throws ErrorException swap!(model, tf, C)

    # # if inserting a node fails, the original nodes should still remain connected
    # @test has_edge(model, bond)
    # @test_throws ErrorException insert_node!(model, bond, Component(:I))
    # @test has_edge(model, bond)
end

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

@testset "Simplifying Junctions" begin
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

@testset "BondGraphNodes" begin
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

@testset "Conversion to Other Graphs" begin
    bg = RCI()
    g = SimpleGraph(bg)
    dg = SimpleDiGraph(bg)

    bg_adj, g_adj, dg_adj = adjacency_matrix.([bg, g, dg])
    @test dg_adj == bg_adj
    @test g_adj == bg_adj + bg_adj' # A + A' forms undirected graph adj matrix
end
