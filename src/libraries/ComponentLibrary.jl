module ComponentLibrary

using ..BondGraphs

using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

const Effort = ModelingToolkitBase.Equality

include("electrical.jl")
include("mechanical.jl")
include("chemical.jl")
include("junctions.jl")

end
