using TestItems
using TestItemRunner

@run_package_tests verbose=true

@testsnippet Setup begin
    using Graphs
    using ModelingToolkit
    using DifferentialEquations: Rosenbrock23
    using Catalyst
    using RecipesBase
    using UnPack

    t = ModelingToolkit.t_nounits
    D = ModelingToolkit.D_nounits

    function RLC()
        r = Component(:R)
        l = Component(:I)
        c = Component(:C)
        kvl = EqualEffort(name=:kvl)

        bg = BondGraph()
        add_node!(bg, [c, l, kvl, r])

        connect!(bg, r, kvl)
        connect!(bg, l, kvl)
        connect!(bg, c, kvl)
        return bg
    end

    function RCI(name=:RCI)
        model = BondGraph(name)
        C = Component(:C)
        R = Component(:R)
        I = Component(:I)
        SS = Component(:SS)
        zero_law = EqualEffort()

        add_node!(model, [C, R, I, SS, zero_law])
        connect!(model, R, zero_law)
        connect!(model, C, zero_law)
        connect!(model, zero_law, I)
        connect!(model, zero_law, SS)

        model
    end

    # cannot use standard notation "var in array" for MTK vars
    var_in(var, dict) = any(iszero.(var .- keys(dict)))

    # sort equations in lex order to make testing equations easier
    sorted_eqs(sys) = sort(equations(sys), by=string)

    function find_subsys(sys, s)
        subsys = ModelingToolkit.get_systems(sys)
        return filter(x -> nameof(x) == s, subsys)[1]
    end
end
