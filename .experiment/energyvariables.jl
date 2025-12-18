import Base: show
using Symbolics

@variables t
D = Differential(t)

# these are symbolic representations and can not yet be "solved"
struct EnergyVariables
    e::Num
    f::Num
    p::Num
    q::Num
end
function EnergyVariables(
        e::Symbol = :e, f::Symbol = :f, p::Symbol = :p, q::Symbol = :q; subscript = nothing)
    if !isnothing(subscript)
        e = Symbol("$(e)_$subscript")
        f = Symbol("$(f)_$subscript")
        p = Symbol("$(p)_$subscript")
        q = Symbol("$(q)_$subscript")
    end
    e, f, p, q = @variables $e(t), $f(t), $p(t), $q(t)
    EnergyVariables(e, f, p, q)
end

effort(ev::EnergyVariables) = ev.e
flow(ev::EnergyVariables) = ev.f
momentum(ev::EnergyVariables) = ev.p
position(ev::EnergyVariables) = ev.q

show(io::IO, ev::EnergyVariables) = print(io, "<$(ev.e),$(ev.f)>")

equations(ev::EnergyVariables) = [D(ev.p) ~ ev.e, D(ev.q) ~ ev.f]

###
ev = EnergyVariables(:u, :v)

equations(ev)

Symbolics.derivative(ev::EnergyVariables) = EnergyVariables(ev.e, ev.f, D(ev.e), D(ev.f))
power(ev::EnergyVariables) = effort(e) * flow(e)
power(e)
Symbolics.derivative(e)
