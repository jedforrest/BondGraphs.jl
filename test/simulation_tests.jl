@variables e(t) f(t) q(t)
@parameters C L R

@testset "Setting variables" begin
    @parameters C=2
    @named ccomp = StaticStorageElement([D(q) ~ f, q ~ C * e], [e], [f], [q]; q=1)
    @named icomp = inductor(L=5)
    @named re = reaction(; normalised=false)
    re.κ = 4

    # getting
    @test ccomp.C == 2 && ccomp.q == 1
    @test icomp.L == 5
    @test re.κ == 4

    # incompatible variables fail
    @test_throws ErrorException re.s = 1
    @test_throws ErrorException ccomp.p = 2
end

@testset "Setting non-numeric control variables" begin
    battery = voltagesource(; name=:se)
    battery.E(t) = sin(2t)
    @test battery.E(0) ≈ 0.0
end

@testset "Simulate RC circuit" begin
    r = resistor(R = 2, name = :R)
    c = capacitor(C = 1, q0=10, name = :C)
    bg = BondGraph([r, c], name=:RC)
    connect!(bg, r, c)

    # true solution
    f_sol(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sys = system(bg)

    sol = solve(ODEProblem(sys, [], tspan))
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[], f_sol(t, 10, 2), atol = 1e-5)
    end

    sol = solve(ODEProblem(sys, [sys.C.q => 5.], tspan))
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f_sol(t, 5, 2), atol = 1e-5)
    end

    sol = solve(ODEProblem(sys, [], tspan, [1., 3.]))
    for t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(t)[1], f_sol(t, 10, 3), atol = 1e-5)
    end
end

@testset "Equivalent resistance (DAE)" begin
    R1 = 1.0
    R2 = 2.0
    Req = R1 + R2
    C = 3.0
    τ = Req * C

    r1 = resistor(name = :r1, R = R1)
    r2 = resistor(name = :r2, R = R2)
    c = capacitor(name = :c; C = C, q0=10)
    kvl = KVL(name = :kvl)
    bg = BondGraph([c, r1, r2, kvl], name=:RRC)
    connect!(bg, c, kvl)
    connect!(bg, kvl, r1)
    connect!(bg, kvl, r2)

    f_sol(x, a, τ) = a * exp(-x / τ)

    tspan = (0.0, 10.0)
    sys = system(bg)

    cq = sys.c.q
    ri = sys.r1.I
    # guesses needed to resolve DAE cyclic conditions
    prob = ODEProblem(sys, [], tspan, guesses=[ri => 1.])
    sol = solve(prob)

    for _t in [0.5, 1.0, 5.0, 10.0]
        @test isapprox(sol(_t)[1], f_sol(_t, 10, τ), atol = 1e-5)
    end
end

@testset "π-filter" begin
    Se = voltagesource(name=:Pin; E = 1)

    Pa = EqualEffort(name = :Pa)
    fa = EqualFlow(name = :fa)
    ca = capacitor(name = :Ca, C = 1, q0 = 1)
    rpa = resistor(name = :Rpa, R = 1)

    Pb = EqualEffort(name = :Pb)
    fb = EqualFlow(name = :fb)
    cb = capacitor(name = :Cb, C = 1, q0 = 2)
    rpb = resistor(name = :Rpb, R = 1)

    fs = EqualFlow(name = :fs)
    l = inductor(name = :L, L = 1, λ0 = 1)
    r = resistor(name = :Rs, R = 1)

    rl = resistor(name = :RL, R = 1)

    bg = BondGraph(name=:π_filter)
    for comp in [Se, Pa, fa, ca, rpa, Pb, fb, cb, rpb, fs, l, r, rl]
        add_comp!(bg, comp)
    end
    connect!(bg, Se, Pa)
    connect!(bg, Pa, fa)
    connect!(bg, fa, ca)
    connect!(bg, fa, rpa)
    connect!(bg, Pa, fs)
    connect!(bg, fs, l)
    connect!(bg, fs, r)
    connect!(bg, fs, Pb)
    connect!(bg, Pb, fb)
    connect!(bg, fb, cb)
    connect!(bg, fb, rpb)
    connect!(bg, Pb, rl)

    tspan = (0, 100.0)
    sys = system(bg)

    constitutive_relations(bg; sub_defaults=true)

    prob = ODEProblem(sys, [], tspan, guesses=[sys.RL.V => 1.5])  # Model is a DAE
    sol = solve(prob, Rosenbrock23())

    (p, Ca, Cb) = (sys.L.λ, sys.Ca.q, sys.Cb.q)
    @test (sol[Ca, 1] ≈ 1) && (sol[Cb, 1] ≈ 2) && (sol[p, 1] == 1)
    @test isapprox(sol[Ca, end], 1.0, atol = 1e-5)
    @test isapprox(sol[Cb, end], 0.5, atol = 1e-5)
    @test isapprox(sol[p, end], 0.5, atol = 1e-5)
end

# TODO Bond Graph Nodes
# @testset "Simulate modular BG" begin
#     r = Component(:R; R = 1)
#     l = Component(:I; L = 1, p = 1)
#     c = Component(:C; C = 1, q = 1)
#     kvl = EqualEffort(name = :kvl)
#     SS1 = SourceSensor(name = :SS1)
#     SS2 = SourceSensor(name = :SS2)

#     bg1 = BondGraph(:RC)
#     add_node!(bg1, [r, c, kvl, SS1])
#     connect!(bg1, r, kvl)
#     connect!(bg1, c, kvl)
#     connect!(bg1, SS1, kvl)
#     bgn1 = BondGraphNode(bg1)

#     bg2 = BondGraph(:L)
#     add_node!(bg2, [l, SS2])
#     connect!(bg2, l, SS2)
#     bgn2 = BondGraphNode(bg2)

#     bg = BondGraph()
#     add_node!(bg, [bgn1, bgn2])
#     connect!(bg, bgn1, bgn2)

#     tspan = (0, 10.0)
#     sol = simulate(bg, tspan)

#     τ = 2
#     ω = sqrt(3) / 2
#     f(t, τ, ω) = exp(-t / τ) * [
#         cos(ω * t) - sqrt(3) * sin(ω * t),
#         cos(ω * t) + sqrt(3) * sin(ω * t)
#     ]

#     for t in [0.0, 0.5, 1.0, 5.0, 10.0]
#         # sort! so that the order of the output is consistent
#         @test isapprox(sort!(sol(t)), sort!(f(t, τ, ω)), atol = 1e-5)
#     end
# end

@testset "Driven Filter Circuit" begin
    model = begin
        @named C = capacitor(C = 1)
        @named R = resistor(R = 1)
        @named Sf = currentsource()
        @named zero_law = KCL()
        model = BondGraph([C, R, zero_law, Sf], name=:RC)
        connect!(model, zero_law, C)
        connect!(model, zero_law, R)
        connect!(model, Sf, zero_law)
        model
    end

    # Simulation parameters
    tspan = (0.0, 5.0)
    u0 = [1]

    # Case 1: constant forcing funciton
    Sf.Is = t -> 3
    @test Sf.Is(t) == 3
    sys = system(model)
    prob = ODEProblem(sys, u0, tspan)
    sol = solve(prob)
    @test isapprox(sol[end], [2.98651], atol = 1e-5)

    # Case 2: regular forcing function
    h(t::Real) = sin(2t)
    Sf.Is = h
    @test isapprox(Sf.Is(pi), 0., atol = 1e-10)
    sys = system(model)
    prob = ODEProblem(sys, u0, tspan)
    sol = solve(prob)
    @test isapprox(sol[end], [0.23625], atol = 1e-5)
end

@testset "Simple Biochemical Simulation" begin
    rn_abc = @reaction_network ABC begin
        1, A + B --> C
    end
    bg_abc = BondGraph(rn_abc)

    constitutive_relations(bg_abc; sub_defaults=true)

    sys = system(bg_abc)
    (A, B, C) = (sys.A.x, sys.B.x, sys.C.x)

    tspan = (0.0, 3.0)
    u0 = [A=>1, B=>2, C=>3]
    prob = ODEProblem(sys, u0, tspan)
    sol = solve(prob)
    @test isapprox(sol[A, end], 1.23606, atol = 1e-5)
    @test isapprox(sol[B, end], 2.23606, atol = 1e-5)
    @test isapprox(sol[C, end], 2.76393, atol = 1e-5)
end

@testset "Stoichiometry Simulation" begin
    rn = @reaction_network A2B begin
        (1, 1), A <--> 2B
    end
    bg = BondGraph(rn)
    sys = system(bg)
    (A, B) = (sys.A.x, sys.B.x)

    prob = ODEProblem(sys, [A=>1.0, B=>0.0], (0.0, 1.0))
    sol = solve(prob)

    # verified by simulation of rn directly
    @test isapprox(sol[A, end], 0.61969, atol = 1e-5)
    @test isapprox(sol[B, end], 0.76062, atol = 1e-5)
end

@testset "Reversible Michaelis-Menten" begin
    rn_mm = @reaction_network MM_reversible begin
        (1, 1), E + S <--> C
        (1, 1), C <--> E + P
    end
    # TODO simplify=true alters the outputted CR
    bg_mm = BondGraph(rn_mm; chemostats = [:S, :P], simplify=false)

    constitutive_relations(bg_mm; sub_defaults=true)

    # setting a new chemostat
    bg_mm[:S].Xs = t -> 2
    @test bg_mm[:S].Xs(t) == 2

    sys = system(bg_mm)
    (E, C) = (sys.E.x, sys.C.x)

    prob = ODEProblem(sys, [E=>1, C=>2], (0.0, 3.0))
    sol = solve(prob)
    @test isapprox(sol[E, end], 1.2, atol = 1e-5)
    @test isapprox(sol[C, end], 1.8, atol = 1e-5)
end

# TODO CONTINUE FROM HERE
# something is wrong with the simplification system
# need to do something else instead of MTK.connect?
# @testset "SERCA (stiff equations)" begin
#     rn_serca = @reaction_network SERCA begin
#         (1, 1), P1 + MgATP <--> P2
#         (1, 1), P2 + H <--> P2a
#         (1, 1), P2 + 2Cai <--> P4
#         (1, 1), P4 <--> P5 + 2H
#         (1, 1), P5 <--> P6 + MgADP
#         (1, 1), P6 <--> P8 + 2Casr
#         (1, 1), P8 + 2H <--> P9
#         (1, 1), P9 <--> P10 + H
#         (1, 1), P10 <--> P1 + Pi
#     end
#     chemostats = [:MgATP, :MgADP, :Pi, :H, :Cai, :Casr]
#     bg_serca = BondGraph(rn_serca; chemostats, simplify=false)

#     reaction_rates = [
#         :Re1 => 0.00053004,
#         :Re2 => 8326784.0537,
#         :Re3 => 1567.7476,
#         :Re4 => 1567.7476,
#         :Re5 => 3063.4006,
#         :Re6 => 130852.3839,
#         :Re7 => 11612934.8748,
#         :Re8 => 11612934.8748,
#         :Re9 => 0.049926
#     ]
#     for (reaction, rate) in reaction_rates
#         bg_serca[reaction].κ = rate
#     end

#     species_affinities = [
#         :P1 => 5263.6085,
#         :P2 => 3803.6518,
#         :P2a => 3110.4445,
#         :P4 => 16520516.1239,
#         :P5 => 0.82914,
#         :P6 => 993148.433,
#         :P8 => 37.7379,
#         :P9 => 2230.2717,
#         :P10 => 410.6048,
#         :Cai => 1.9058,
#         :Casr => 31.764,
#         :MgATP => 244.3021,
#         :MgADP => 5.8126e-7,
#         :Pi => 0.014921,
#         :H => 1862.5406
#     ]
#     for (species, affinity) in species_affinities
#         bg_serca[species].K = affinity
#     end

#     chemostat_amounts = [
#         :Cai => t -> 0.0057,
#         :Casr => t -> (0.05 + 0.01t)*2.28,
#         :H => t -> 0.004028,
#         :MgADP => t -> 1.3794,
#         :MgATP => t -> 3.8,
#         :Pi => t -> 570
#     ]
#     for (chemostat, amount) in chemostat_amounts
#         bg_serca[chemostat].Xs = amount
#     end

#     initial_conditions = [
#         :P1 => 0.000483061870385487,
#         :P2 => 0.0574915174273067,
#         :P2a => 0.527445119834607,
#         :P4 => 1.51818391164022e-09,
#         :P5 => 0.000521923287622898,
#         :P6 => 7.80721128535043e-05,
#         :P8 => 0.156693953834181,
#         :P9 => 0.149232225342376,
#         :P10 => 0.108044124948978
#     ]
#     # for (species, ic) in initial_conditions
#     #     getproperty(bg_serca, species).q = ic
#     # end

#     tspan = (0.0, 200.0)
#     sys = system(bg_serca)
#     prob = ODEProblem(sys, initial_conditions, tspan)
#     sol = simulate(bg_serca, tspan; solver = Rosenbrock23());

#     # calculated using the same model, verified by plot from BGT tutorial
#     sys = ODESystem(bg_serca, simplify_eqs = false)
#     real_solution = Dict(
#         sys.P1.x => 4.4404656222265794e-5,
#         sys.P2.x => 0.09777422826977565,
#         sys.P2a.x => 0.8970112784324162,
#         sys.P4.x => 2.6596475539704174e-9,
#         sys.P5.x => 0.0009426424413096248,
#         sys.P6.x => 0.001015195974904865,
#         sys.P8.x => 0.001212098675876874,
#         sys.P9.x => 0.0011543788312157496,
#         sys.P10.x => 0.0008357702283899367
#     )

#     for (var, real_sol) in real_solution
#         @test isapprox(sol[var, end], real_sol, atol = 1e-5)
#     end
# end
