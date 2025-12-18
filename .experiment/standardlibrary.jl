module Library

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics: scalarize

const Effort = ModelingToolkit.Equality

############################################################

@connector PowerPort begin
    @structural_parameters begin
        N = 1
    end
    @variables begin
        e(t)[1:N] = 0.0, [connect = Effort]
        f(t)[1:N] = 0.0, [connect = Flow]
    end
end

############################################################

@connector PowerVars begin
    @structural_parameters begin
        N = 1
    end
    @variables begin
        e(t)[1:N] = 0.0, [connect = Effort]
        f(t)[1:N] = 0.0, [connect = Flow]
    end
end
@named pv = PowerVars(N = 1)

@mtkmodel StaticStorage begin
    @description "Static Storage Element"
    @structural_parameters begin
        phi = (e, q, C) -> q ~ C * e
    end
    @extend PowerVars()
    @parameters begin
        C
    end
    @variables begin
        q(t)[1:N] = 0.0
    end
    @equations begin
        D(q) ~ f
        phi(e, q, C)
    end
end

############################################################
# Electrical

@mtkmodel ElectricalVariables begin
    @variables begin
        V(t), [connect = Effort]
        I(t), [connect = Flow]
    end
end

@mtkmodel Capacitor begin
    @description "Generalised Linear Capacitor"
    @extend ElectricalVariables()
    @variables begin
        q(t)
    end
    @defaults begin
        q => 0.0
    end
    @parameters begin
        C = 1
    end
    @equations begin
        D(q) ~ I
        q ~ C * V
    end
end

@mtkmodel Resistor begin
    @description "Generalised Linear Resistor"
    @extend ElectricalVariables()
    @parameters begin
        R = 1
    end
    @equations begin
        V ~ R * I
    end
end

@mtkmodel Inductor begin
    @description "Generalised Linear Inductor"
    @extend ElectricalVariables()
    @variables begin
        p(t)
    end
    @defaults begin
        p => 0.0
    end
    @parameters begin
        L = 1
    end
    @equations begin
        D(p) ~ V
        p ~ L * I
    end
end

@mtkmodel EffortSource begin
    @description "Effort Source"
    @extend ElectricalVariables()
    @parameters begin
        E = 1.0
    end
    @equations begin
        V ~ E
    end
end

@mtkmodel FlowSource begin
    @description "Flow Source"
    @extend ElectricalVariables()
    @parameters begin
        F = 1.0
    end
    @equations begin
        I ~ F
    end
end

@named cap = Capacitor(C = 1.0)
@named res = Resistor(R = 2.0)
@named ind = Inductor(L = 3.0)

@named esrc = EffortSource(E = 5.0)
@named fsrc = FlowSource(F = 2.0)

############################################################
# Junctions
############################################################

@mtkmodel Transformer begin
    @description "Ideal Transformer"
    @components begin
        port = [ElectricalVariables() for i in 1:2]
    end
    @parameters begin
        n = 1
    end
    @equations begin
        port[1].V ~ n * port[2].V
        port[2].I ~ n * port[1].I
    end
end

@mtkmodel Gyrator begin
    @description "Ideal Gyrator"
    @components begin
        port = [ElectricalVariables() for i in 1:2]
    end
    @parameters begin
        r = 1
    end
    @equations begin
        port[1].V ~ r * port[2].I
        port[2].V ~ r * port[1].I
    end
end

############################################################

check_num_ports(N) = N >= 2 || error("Junction must have at least 2 ports ($N ports given)")

@mtkmodel ZeroJunction begin
    @description "0-Junction (EqualEffort)"
    begin
        check_num_ports(N)
    end
    @extend PowerPort(; N)
    @equations begin
        scalarize(sum(f)) ~ 0
        scalarize([e[1] ~ e_i for e_i in e[2:end]])
    end
end

@mtkmodel OneJunction begin
    @description "1-Junction (EqualFlow)"
    begin
        check_num_ports(N)
    end
    @extend PowerPort(; N)
    @equations begin
        scalarize(sum(e)) ~ 0
        scalarize([f[1] ~ f_i for f_i in f[2:end]])
    end
end

@named tfmr = Transformer(n = 3.0)
@named gyr = Gyrator(r = 4.0)

@named zjunc = ZeroJunction(N = 4)
@named ojunc = OneJunction(N = 10)
equations(zjunc)
equations(ojunc)

############################################################################################
# Biochemical

const _R = 8.314
const _T = 310.0

@mtkmodel ChemicalSpecies begin
    @description "Chemical Species"
    @components begin
        port = PowerPort()
    end
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
        port.e[1] ~ R * T * log(K * x)    # chemical potential
        port.f[1] ~ D(x)          # flow is time derivative of amount
    end
end

@mtkmodel Reaction begin
    @description "Chemical Reaction (Re)"
    @components begin
        reactants = PowerPort()
        products = PowerPort()
    end
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
        Af ~ exp(reactants.e[1] / (R * T))
        Ar ~ exp(products.e[1] / (R * T))
        reactants.f[1] ~ products.f[1]
        reactants.f[1] ~ r * (Af - Ar)
    end
end

@mtkmodel ChemicalSource begin
    @description "Chemical Source (Se)"
    @components begin
        port = PowerPort()
    end
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
        port.e[1] ~ R * T * log(K * X)
    end
end

@named ce = ChemicalSpecies(K = 2.0)
@named re = Reaction(r = 0.5)
@named se = ChemicalSource(X = 1.0)

end
