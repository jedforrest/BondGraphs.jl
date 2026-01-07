### Electrical units
V, I, λ, q = power_variables(e="V", f="I", p="λ", q="q")

# Resistor
@parameters R=1
resistor = DissipatorElement([V ~ R * I], [V], [I])

# Capacitor
@parameters C=1
capacitor = StaticStorageElement([D(q) ~ I, q ~ C * V], [V], [I], [q])

# Inductor
@parameters L=1
inductor = DynamicStorageElement([D(λ) ~ V, λ ~ L * I], [V], [I], [λ])

# Sources
@parameters E=1
voltagesource = EffortSource([V ~ E], [V], [I])

@parameters F=1
currentsource = FlowSource([F ~ I], [V], [I])

# Junctions
zerojunction = EqualEffort()
onejunction = EqualFlow()
