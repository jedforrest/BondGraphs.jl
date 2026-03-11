const _R = 1
const _T = 1
@parameters R=_R T=_T
R = GlobalScope(R)
T = GlobalScope(T)

# Chemical Species
function chemicalspecies(; name)
    @parameters K=1
    u, v, x = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
        x(t)
    end
    StaticStorageElement([D(x) ~ v, u ~ R*T*log(K*x)], [u], [v], [x]; name)
end

# Reaction
function reaction(; name)
    @parameters κ=1
    u, v = @variables begin
        μ(t)[1:2], [connect = Effort]
        ν(t)[1:2], [connect = Flow]
    end
    DissipatorElement([v[1] ~ - v[2], v[1] ~ κ * (exp(u[1] / (R*T)) - exp(u[2] / (R*T)))], u, v; name)
end

# Chemical Source (chemostat)
function chemostat(; name)
    @parameters K=1 X=1
    u, v = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
    end
    EffortSource([u ~ R*T*log(K*X)], [u], [v]; name)
end

# Stoichiometry
function stoichiometry(n; name=Symbol("TF{$n}"))
    Transformer(n; name)
end
