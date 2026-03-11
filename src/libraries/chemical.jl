@constants begin
    T=1
    R=1
end
const RT = R * T

# Chemical Species
function chemicalspecies(; name)
    @parameters K=1
    u, v, x = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
        x(t)
    end
    StaticStorageElement([D(x) ~ v, u ~ RT*log(K*x)], [u], [v], [x]; name)
end

# Reaction
function reaction(; name)
    @parameters κ=1
    u, v = @variables begin
        μ(t)[1:2], [connect = Effort]
        ν(t)[1:2], [connect = Flow]
    end
    DissipatorElement([v[1] ~ - v[2], v[1] ~ κ * (exp(u[1] / RT) - exp(u[2] / RT))], u, v; name)
end

# Chemical Source (chemostat)
function chemostat(; name)
    @parameters K=1 X=1
        u, v = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
    end
    EffortSource([u ~ RT*log(K*X)], [u], [v]; name)
end

# Stoichiometry
function stoichiometry(n; name=Symbol("TF{$n}"))
    Transformer(n; name)
end
