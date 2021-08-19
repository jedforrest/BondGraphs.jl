using BondGraphs
using Test

c1 = Component(:C, "C1")
c2 = Component(:C, "C2")
c1 == c2

fieldnames(Bond)

b1 = Bond(c1,c2)
Bond

import LightGraphs; const lg = LightGraphs
lg.dst(b1)

supertype(Component)

@testset "BondGraphs.jl" begin
    # Write your tests here.
end
