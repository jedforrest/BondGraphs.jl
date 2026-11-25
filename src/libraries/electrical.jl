### Electrical units
@variables begin
    V(t), [connect = Effort]
    I(t), [connect = Flow]
    λ(t),
    q(t)
end

# Resistor
function resistor(; R=1, name)
    @parameters R=R
    Component(Resistance(), [V ~ R * I], [V], [I]; name)
end

# Capacitor
function capacitor(; C=1, q=0, name)
    @parameters C=C
    Component(Compliance(), [D(q) ~ I, q ~ C * V], [V], [I], [q]; name, q=q)
end

# Inductor
function inductor(; L=1, λ=0, name)
    @parameters L=L
    Component(Inertance(), [D(λ) ~ V, λ ~ L * I], [V], [I], [λ]; name, λ=λ)
end

# Sources
function voltagesource(; E=t->1, name)
    @parameters Vs(::Real) = E
    Component(EffortSource(), [V ~ Vs(t)], [V], [I]; name)
end

function currentsource(; F=t->1, name)
    @parameters Is(::Real) = F
    Component(FlowSource(), [I ~ Is(t)], [V], [I]; name)
end
