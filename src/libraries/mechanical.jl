### Mechanical units
@variables begin
    F(t), [connect = Effort]
    v(t), [connect = Flow]
    p(t),
    x(t)
end

# Damper
function damper(; b=1, name)
    @parameters b=b
    Component(Resistance(), [F ~ b * v], [F], [v]; name)
end

# Spring
function spring(; k=1, x=0, name)
    @parameters k=k
    Component(Compliance(), [D(x) ~ F, x ~ k * F], [F], [x]; name, x=x)
end

# Inertial Mass
function mass(; m=1, p=0, name)
    @parameters m=m
    Component(Inertance(), [D(p) ~ F, p ~ m * F], [F], [p]; name, p=p)
end

# Sources
function externalforce(; F=t->1, name)
    @parameters Fs(::Real) = F
    Component(EffortSource(), [F ~ Fs(t)], [F], [v]; name)
end

function externalvelocity(; v=t->1, name)
    @parameters vs(::Real) = v
    Component(FlowSource(), [v ~ vs(t)], [F], [v]; name)
end
