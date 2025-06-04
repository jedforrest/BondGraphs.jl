using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

@mtkbuild fol = FOL()
prob = ODEProblem(fol, [], (0.0, 10.0), [])
sol = solve(prob)

plot(sol)

FOL
supertypes(typeof(FOL))


@variables x(t)=1 y(t)=0 vx(t)=0 vy(t)=2

continuous_events = [[x ~ 0] => [vx ~ -vx]
                     [y ~ -1.5, y ~ 1.5] => [vy ~ -vy]]

@mtkbuild ball = ODESystem(
    [
        D(x) ~ vx,
        D(y) ~ vy,
        D(vx) ~ -9.8 - 0.1vx, # gravity + some small air resistance
        D(vy) ~ -0.1vy
    ], t; continuous_events)


@connector function HydraulicPort(; p_int, name)
    pars = @parameters begin
        ρ
        β
        μ
    end

    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end

    ODESystem(Equation[], t, vars, pars; name, defaults = [dm => 0])
end

@named h_port = HydraulicPort(p_int = 1)
fieldnames(typeof(h_port))
h_port.connector_type
