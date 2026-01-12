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
    j0 = EqualEffort(name=:𝟎)
    j1 = EqualFlow(name=:𝟏)
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

@testset "Inserting Nodes" begin
    using BondGraphs: transformer, KCL
    bg = RCI()

    bond_c_kvl = bg[:c, :kvl]
    bond_r_kvl = bg[:r, :kvl]

    tf = transformer(2)
    insert_comp!(bg, bond_c_kvl, tf)

    @named kcl = KCL()
    insert_comp!(bg, bond_r_kvl, kcl)

    @test tf in components(bg) && kcl in components(bg)
    @test nv(bg) == 7
    @test ne(bg) == 6
end

@testset "Construction Failure" begin
    using BondGraphs: resistor, capacitor, inductor, KVL
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
    tf = transformer(2)
    add_comp!(bg, tf)
    @test_logs (:warn, "New comp must have a greater or equal number of ports than the old comp") swap!(bg, tf, ccomp)

    # if inserting a node fails, the original nodes should still remain connected
    @named icomp = inductor()
    @test_throws ErrorException insert_comp!(bg, b1, icomp)
    @test b1 in bg
end

### TODO CONTINUE FROM HERE
@testset "Simplifying Junctions" begin
    bg = RCI()
    C, R, I, SS, J1 = bg[:c], bg[:r], bg[:i], bg[:v], bg[:kvl]

    # adding redundancies
    J1_new_1 = EqualFlow(; name = :new1_1)
    J1_new_2 = EqualFlow(; name = :new1_2)
    insert_comp!(bg, (C, J1), J1_new_1)
    insert_comp!(bg, (R, J1), J1_new_2)
    connect!(bg, J1_new_1, J1_new_2)

    J0_new_1 = EqualEffort(; name = :new0_1)
    J0_new_2 = EqualEffort(; name = :new0_2)
    add_comp!(bg, J0_new_1)
    connect!(bg, J1, J0_new_1)
    insert_comp!(bg, (J1, SS), J0_new_2)

    # Removing junction redundancies
    juncs = junctions(bg)
    zero_juncs = filter(j -> j isa EqualEffort, juncs)
    one_juncs = filter(j -> j isa EqualFlow, juncs)
    @test length(zero_juncs) == 2 && length(one_juncs) == 3

    simplify_junctions!(bg, squash_identical = false)
    zero_juncs = filter(j -> j isa EqualEffort, junctions(bg))
    @test length(zero_juncs) == 0
    @test nv(bg) == 7
    @test ne(bg) == 7

    # Squashing junction duplicates into a single junction
    # FIXME doesn't like swapping with another component already in the graph
    # need to write code instead of using 'swap!'
    # simplify_junctions!(bg)
    # one_juncs = filter(j -> j isa EqualFlow, junctions(bg))
    # @test length(one_juncs) == 1
    # @test nv(bg) == 5
    # @test ne(bg) == 4
end
