### Electrical units
V, I, λ, q = power_variables(e="V", f="I", p="λ", q="q")

# Resistor
function resistor(; name)
    @parameters R=1
    DissipatorElement([V ~ R * I], [V], [I]; name)
end

# Capacitor
function capacitor(; name)
    @parameters C=1
    StaticStorageElement([D(q) ~ I, q ~ C * V], [V], [I], [q]; name)
end

# Inductor
function inductor(; name)
    @parameters L=1
    DynamicStorageElement([D(λ) ~ V, λ ~ L * I], [V], [I], [λ]; name)
end

# Sources
function voltagesource(; name)
    @parameters E=1
    EffortSource([V ~ E], [V], [I]; name)
end

function currentsource(; name)
    @parameters F=1
    FlowSource([F ~ I], [V], [I]; name)
end

KCL(; name) = EqualEffort(; name)
KVL(; name) = EqualFlow(; name)
