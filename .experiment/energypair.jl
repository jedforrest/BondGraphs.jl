import Base: show
using Symbolics

# these are symbolic representations and can not yet be "solved"
struct EnergyPair
    e::Num
    f::Num
end
function EnergyPair(e::Symbol=:e, f::Symbol=:f)
    _, e, f = @variables t, $e(t), $f(t)
    EnergyPair(e, f)
end

effort(ep::EnergyPair) = ep.e
flow(ep::EnergyPair) = ep.f

show(io::IO, ep::EnergyPair) = print(io, "<$(effort(ep)),$(flow(ep))>")

power(ep::EnergyPair) = effort(ep) * flow(ep)

ep = EnergyPair(:u, :v)
