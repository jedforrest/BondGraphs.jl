### Electrical units
V, I, λ, q = power_variables(e="V", f="I", p="λ", q="q")

# Resistor
function resistor(; R=1, name)
    @parameters R=R
    DissipatorElement([V ~ R * I], [V], [I]; name)
end

# Capacitor
function capacitor(; C=1, q0=0, name)
    @parameters C=C
    StaticStorageElement([D(q) ~ I, q ~ C * V], [V], [I], [q]; name, q=q0)
end

# Inductor
function inductor(; L=1, λ0=0, name)
    @parameters L=L
    DynamicStorageElement([D(λ) ~ V, λ ~ L * I], [V], [I], [λ]; name, λ=λ0)
end

# Sources
function voltagesource(; E=1, name)
    @parameters E(t)=E
    EffortSource([V ~ E], [V], [I]; name)
end

function currentsource(; F=1, name)
    @parameters F(t)=F
    FlowSource([I ~ F], [V], [I]; name)
end

# Transformers
function transformer(n; name=Symbol("TF{$n}"))
    Transformer(n; name)
end

function gyrator(r; name=Symbol("GY{$n}"))
    Gyrator(r; name)
end

# Junctions
KCL(; name) = EqualEffort(; name)
KVL(; name) = EqualFlow(; name)
