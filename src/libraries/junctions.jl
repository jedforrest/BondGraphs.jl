# Transformers
function transformer(n; name=Symbol("TF{$n}"))
    @variables begin
        e(t)[1:2], [connect = Effort]
        f(t)[1:2], [connect = Flow]
    end
    eqs = [
        e[1] ~ n * e[2],
        f[2] ~ -n * f[1]
    ]
    Component(Transformer(), eqs, [e[1], e[2]], [f[1], f[2]]; name)
end

function gyrator(r; name=Symbol("GY{$r}"))
    @variables begin
        e(t)[1:2], [connect = Effort]
        f(t)[1:2], [connect = Flow]
    end
    eqs = [
        e[1] ~ r * f[2],
        e[2] ~ -r * f[1]
    ]
    Component(Gyrator(), eqs, [e[1], e[2]], [f[1], f[2]]; name)
end

# Junctions
zerojunction(; name) = Component(ZeroJunction(); name)
onejunction(; name) = Component(OneJunction(); name)
