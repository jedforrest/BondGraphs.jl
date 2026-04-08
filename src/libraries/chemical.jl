# NOTE: Chemical reactions are nonlinear, so bond direction matters
# For multiple species, the bonds into the junction(s) on each side of the reaction
# must have the same direction (inbound or outbound)
# e.g. A + B <=> C
# A -> 1
# B -> 1
# 1 -> Re
# Re -> C

const _R = 8.314
const _T = 310
@parameters R=_R T=_T
R = GlobalScope(R)
T = GlobalScope(T)

# Chemical Species
function chemicalspecies(; K=1, normalised=true, name)
    @parameters K=K
    u, v, x = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
        x(t)
    end
    eqns = if normalised
        [D(x) ~ v, u ~ log(K*x)]
    else
        [D(x) ~ v, u ~ R*T*log(K*x)]
    end
    StaticStorageElement(eqns, [u], [v], [x]; name)
end

# Reaction
function reaction(; κ=1, normalised=true, name)
    @parameters κ=κ
    u, v = @variables begin
        μ(t)[1:2], [connect = Effort]
        ν(t)[1:2], [connect = Flow]
    end
    eqns = if normalised
        [v[1] ~ - v[2], v[1] ~ κ * (exp(u[1]) - exp(u[2]))]
    else
        [v[1] ~ - v[2], v[1] ~ κ * (exp(u[1] / (R*T)) - exp(u[2] / (R*T)))]
    end
    DissipatorElement(eqns, u, v; name)
end

# Chemical Source (chemostat)
function chemostat(; K=1, X=1, normalised=true, name)
    @parameters K=K X=X
    u, v = @variables begin
        μ(t), [connect = Effort]
        ν(t), [connect = Flow]
    end
    eqns = if normalised
        [D(x) ~ v, u ~ log(K*X)]
    else
        [D(x) ~ v, u ~ R*T*log(K*X)]
    end
    EffortSource(eqns, [u], [v]; name)
end

# Stoichiometry
function stoichiometry(n; name=Symbol("TF{$n}"))
    Transformer(n; name)
end
