import Base: show, iterate
using ModelingToolkit
using Symbolics

# these are symbolic representations and can not yet be "solved"
struct EffortFlowPair
    e::Num
    f::Num
end
function EffortFlowPair(e::Symbol=:e, f::Symbol=:f)
    _, e_var, f_var = @variables t, $e(t), $f(t)
    EffortFlowPair(e_var, f_var)
end

effort(ef::EffortFlowPair) = ef.e
flow(ef::EffortFlowPair) = ef.f

# allows iteration/unpacking over the struct
iterate(ef::EffortFlowPair) = (ef.e, true)
iterate(ef::EffortFlowPair, state) = state ? (ef.f, false) : nothing


# # TODO define p and q using symbolics integral
# integrate(var::Num) = Integral(Symbolics.VarDomainPairing(var, [0, t]))
# momentum(ef::EnergyPair) = integrate(effort(ep))
# position(ef::EnergyPair) = integrate(flow(ep))

power(ef::EffortFlowPair) = effort(ef) * flow(ef)

show(io::IO, ef::EffortFlowPair) = print(io, "<$(effort(ef)), $(flow(ef))>")

####################################

@connector EffortFlow begin
    e(t)
    f(t), [connect = Flow]
end

####################################
ef = EffortFlowPair()
e_num, f_num = ef
uv = EffortFlowPair(:u, :v)

effort(ef)
flow(ef)
power(ef)

efs = map(i -> EffortFlowPair(Symbol("e$i"), Symbol("f$i")), 1:10)
effort.(efs)
flow.(efs)
power.(efs)
