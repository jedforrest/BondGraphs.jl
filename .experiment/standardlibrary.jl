module Library

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics: scalarize

const Effort = ModelingToolkit.Equality

@connector PowerPort begin
    @structural_parameters begin
        N = 1
    end
    @variables begin
        e(t)[1:N] = 0., [connect = Effort]
        f(t)[1:N] = 0., [connect = Flow]
    end
end

@named ps = PowerPort(N=2)

############################################################

@mtkmodel Capacitor begin
    @description "Generalised Linear Capacitor"
    @extend PowerPort(N=1)
    @parameters begin
        C = 1.0
    end
    @variables begin
        q(t)
    end
    @defaults begin
        q => 0.0
    end
    @equations begin
        D(q) ~ f[1]
        q ~ C * e[1]
    end
end

@mtkmodel Resistor begin
    @description "Generalised Linear Resistor"
    @extend PowerPort(N=1)
    @parameters begin
        R = 1.0
    end
    @equations begin
        e[1] ~ R * f[1]
    end
end

@mtkmodel Inductor begin
    @description "Generalised Linear Inductor"
    @extend PowerPort(N=1)
    @parameters begin
        L = 1.0
    end
    @variables begin
        p(t)
    end
    @defaults begin
        p => 0.0
    end
    @equations begin
        p ~ L * f[1]
        D(p) ~ e[1]
    end
end

@mtkmodel EffortSource begin
    @description "Effort Source"
    @extend PowerPort(N=1)
    @parameters begin
        E = 1.0
    end
    @equations begin
        e[1] ~ E
    end
end

@mtkmodel FlowSource begin
    @description "Flow Source"
    @extend PowerPort(N=1)
    @parameters begin
        F = 1.0
    end
    @equations begin
        f[1] ~ F
    end
end

@named cap = Capacitor(C=1.0)
@named res = Resistor(R=2.0)
@named ind = Inductor(L=3.0)

@named esrc = EffortSource(E=5.0)
@named fsrc = FlowSource(F=2.0)

############################################################
# Junctions

check_num_ports(N) = N >= 2 || error("Junction must have at least 2 ports (N=$N)")

@mtkmodel Transformer begin
    @description "Ideal Transformer"
    @extend PowerPort(N=2)
    @parameters begin
        n = 1.0
    end
    @equations begin
        e[1] ~ n * e[2]
        f[2] ~ n * f[1]
    end
end

@mtkmodel Gyrator begin
    @description "Ideal Gyrator"
    @extend PowerPort(N=2)
    @parameters begin
        r = 1.0
    end
    @equations begin
        e[1] ~ r * f[2]
        e[2] ~ r * f[1]
    end
end

@mtkmodel ZeroJunction begin
    @description "0-Junction (EqualEffort)"
    @structural_parameters begin
        Nports = 2
    end
    begin
        check_num_ports(Nports)
    end
    @extend PowerPort(N=Nports)
    @equations begin
        scalarize(sum(f)) ~ 0
        scalarize([e[1] ~ e_i for e_i in e[2:end]])
    end
end

@mtkmodel OneJunction begin
    @description "1-Junction (EqualFlow)"
    @structural_parameters begin
        Nports = 2
    end
    begin
        check_num_ports(Nports)
    end
    @extend PowerPort(N=Nports)
    @equations begin
        scalarize(sum(e)) ~ 0
        scalarize([f[1] ~ f_i for f_i in f[2:end]])
    end
end


@named tfmr = Transformer(n=3.0)
@named gyr = Gyrator(r=4.0)

@named zjunc = ZeroJunction(Nports=4)
@named ojunc = OneJunction(Nports=10)
equations(zjunc)
equations(ojunc)

############################################################################################
# Biochemical

const _R = 8.314
const _T = 310.0

@mtkmodel ChemicalSpecies begin
    @description "Chemical Species"
    @extend PowerPort(N=1)
    @constants begin
        R = _R
        T = _T
    end
    @parameters begin
        K = 1.0      # thermodynamic constant
    end
    @variables begin
        x(t)         # amount or concentration
    end
    @defaults begin
        x => 0.0
    end
    @equations begin
        e[1] ~ R * T * log(K * x)    # chemical potential
        f[1] ~ D(x)          # flow is time derivative of amount
    end
end

@mtkmodel Reaction begin
    @description "Chemical Reaction (Re)"
    @extend PowerPort(N=2)
    @constants begin
        R = _R
        T = _T
    end
    @parameters begin
        r = 1.0
    end
    @variables begin
        Af(t)  # forward affinity
        Ar(t)  # reverse affinity
    end
    @equations begin
        f[1] ~ f[2]
        Af ~ exp(e[1] / (R * T))
        Ar ~ exp(e[2] / (R * T))
        f[1] ~ r * (Af - Ar)
    end
end

@mtkmodel ChemicalSource begin
    @description "Chemical Source (Se)"
    @extend PowerPort(N=1)
    @constants begin
        R = _R
        T = _T
    end
    @parameters begin
        K = 1.0      # thermodynamic constant
    end
    @parameters begin
        X = 0.0    # chemical potential
    end
    @equations begin
        e[1] ~ R * T * log(K * X)
    end
end

@named ce = ChemicalSpecies(K=2.0)
@named re = Reaction(r=0.5)
@named se = ChemicalSource(X=1.0)

end
