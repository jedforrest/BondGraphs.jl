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
function capacitor(; C=1, q0=0, name)
    @parameters C=C
    Component(Compliance(), [D(q) ~ I, q ~ C * V], [V], [I], [q]; name, q=q0)
end

# Inductor
function inductor(; L=1, λ0=0, name)
    @parameters L=L
    Component(Inertance(), [D(λ) ~ V, λ ~ L * I], [V], [I], [λ]; name, λ=λ0)
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

# Transformers
function transformer(n; name=Symbol("TF{$n}"))
    @variables begin
        v(t)[1:2], [connect = Effort]
        i(t)[1:2], [connect = Flow]
    end
    eqs = [
        v[1] ~ n * v[2],
        i[2] ~ -n * i[1]
    ]
    Component(Transformer(), eqs, [v[1], v[2]], [i[1], i[2]]; name)
end

function gyrator(r; name=Symbol("GY{$r}"))
    @variables begin
        v(t)[1:2], [connect = Effort]
        i(t)[1:2], [connect = Flow]
    end
    eqs = [
        v[1] ~ r * i[2],
        v[2] ~ -r * i[1]
    ]
    Component(Gyrator(), eqs, [v[1], v[2]], [i[1], i[2]]; name)
end

# Junctions
KCL(; name) = Component(ZeroJunction(); name)
KVL(; name) = Component(OneJunction(); name)
