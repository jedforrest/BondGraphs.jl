time_dependent_var(x) = 
    Num(Variable{ModelingToolkit.FnType{Tuple{Any}, Real}}(x))(TIME)
real_var(x) = ModelingToolkit.toparam(Num(Variable{Real}(x)))

internal_effort(idx) = time_dependent_var(Symbol("e_$idx"))
internal_flow(idx) = time_dependent_var(Symbol("f_$idx"))
external_effort(idx) = time_dependent_var(Symbol("E_$idx"))
external_flow(idx) = time_dependent_var(Symbol("F_$idx"))
internal_state(idx) = time_dependent_var(Symbol("x_$idx"))
control_variable(idx) = real_var(Symbol("u_$idx"))

# Equations
equations(n::Component) = n.equations
function equations(m::EqualEffort)
    if m.degree == 0
        return Vector{Equation}([])
    end
    e_vars = external_effort.(1:m.degree)
    f_vars = external_flow.(1:m.degree)

    flow_constraint = [0 ~ sum(f_vars)]
    effort_constraints = [0 ~ e_vars[1] - e for e in e_vars[2:end]]
    return vcat(flow_constraint,effort_constraints)
end
function equations(m::EqualFlow)
    if m.degree == 0
        return Vector{Equation}([])
    end
    e_vars = external_effort.(1:m.degree)
    f_vars = external_flow.(1:m.degree)

    weighted_e = m.weights.*e_vars
    weighted_f = m.weights.*f_vars

    effort_constraint = [0 ~ sum(weighted_e)]
    flow_constraints = [0 ~ weighted_f[1] - f for f in weighted_f[2:end]]
    return vcat(effort_constraint,flow_constraints)
end