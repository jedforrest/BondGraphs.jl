# This will replicate the existing 'getting started' tutorial with the new syntax

include("../ontology.jl")

include("../standardlibrary.jl")
using .Library

@named ppt = PowerPort()

@named res = Library.Resistor()
@named cap = Library.Capacitor()

# R: e <=> f
r = DissipatorElement(res)
# C: e <=> q
c = StaticStorageElement(cap)

System(equations(r), t, name = :r)

@named rcomp = Component(r)
@named ccomp = Component(c)
@named zcomp = Component(EqualEffort())

b1 = Bond(rcomp, zcomp)
b2 = Bond(ccomp, zcomp)

rcomp.ports
zcomp.ports

bg = BondGraph("RC Circuit", [rcomp, ccomp, zcomp], [b1, b2])

system(rcomp)
system(ccomp)
sys0 = system(zcomp)
equations(expand_connections(sys0))

bgsys = system(bg; simplify = false)
hierarchy(bgsys)

sys2 = structural_simplify(bgsys)  # will become mtkcompile in later update

unknowns(sys2)
equations(sys2)
observed(sys2)
full_equations(sys2)
ModelingToolkit.iscomplete(sys2)

#############
prob = ODEProblem(sys2, [sys2.ccomp.q => 5], (0.0, 10.0), [
    sys2.rcomp.R => 2, sys2.ccomp.C => 1])
prob.ps
sol = solve(prob, Tsit5())
plot(sol)

#############
using Graphs
g = graph(bg)
incidence_matrix(g)
graphplot(bg)
#############

@named flow = Library.FlowSource()
A = FlowSource(flow)

@named fcomp = Component(A)
b3 = Bond(fcomp, zcomp)
bg2 = BondGraph("RCA Circuit", [rcomp, ccomp, zcomp, fcomp], [b1, b2, b3])
graphplot(bg2)

bg2_sys = system(bg2; simplify = true)

full_equations(bg2_sys)
