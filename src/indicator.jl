export indicator

"""
    indicator(y, b; indregion=:above, ineq=:nonstrict)

Indicator function of default form `1\{y\\geq b\}`. The inequality can be made strict and the direction of inequality can be reversed.

# Arguments
- `indregion=:above` (default): calculates ``1\{y\\geq b\}`` or ``1\{y>b\}``.
- `indregion=:below`: calculates ``1\{y\\leq b\}`` or ``1\{y<b\}``.
- `ineq=:nonstrict` (default): calculates ``1\{y\\geq b\}`` or ``1\{y\\leq b\}``.
- `ineq=:strict`: calculates ``1\{y>b\}`` or ``1\{y<b\}``.
"""
function indicator(y, b; indregion = :above, ineq=:nonstrict)
    # Indicator function. The default mode ('above') returns 𝟏{y ≥ b}. The other available mode, 'below', returns 𝟏{y ≤ b}. Set 'ineq' to 'strict' to get 𝟏{y > b} or 𝟏{y < b}.
    indregion = Symbol(indregion)
    ineq = Symbol(ineq)
    if ineq == :nonstrict
        if indregion == :above
            y ≥ b ? (return 1) : (return 0)
        end
        if indregion == :below
            y ≤ b ? (return 1) : (return 0)
        end
    end
    if ineq == :strict
        if indregion == :above
            y > b ? (return 1) : (return 0)
        end
        if indregion == :below
            y < b ? (return 1) : (return 0)
        end
    end
    regionlist = [:above, :below]
    if (indregion .== regionlist) == zeros(2)
        error("Unknown indicator region: $indregion")
    end
    ineqlist = [:strict, :nonstrict]
    if (indregion .== ineqlist) == zeros(2)
        error("Unknown inequality type: $ineq")
    end
end
