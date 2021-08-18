using ThresholdStability
using Plots, Measures
using LaTeXStrings
pyplot()

E1, E2, E3, E4 = [1 0.; 0 1.], [1 0.; 0 -1.], [-1 0.; 0 1.], [-1 0.; 0 -1.]
D1 = zeros(1,2); D2, D3, D4 = copy(D1), copy(D1), copy(D1)
X = [[E1, D1], [E2, D2], [E3, D3], [E4, D4]]

function AR2_to_companion(ϕs, ϕ_stars)
    # Returns the set of matrices 𝐀 from putting AR(2) model
    #   yₜ* = ϕ₁*yₜ₋₁* + ϕ₁yₜ₋₁ + ϕ₂*yₜ₋₂* + ϕ₂yₜ₋₂ +
    #   yₜ = max{yₜ*, 0}
    #   in form
#            𝐲ₜ = 𝐀𝐲ₜ₋₁ + 𝛆ₜ
#    where   𝐲ₜ = [yₜ*, yₜ₋₁*, yₜ₋₁]; 𝛆ₜ = [ϵₜ, 0, 0]
#    and     𝐀 = [ϕ₁*+ϕ₁𝐈{yₜ₋₁* ≥ 0}   ϕ₂*  ϕ₂
#                        1             0    0
#                    𝐈{yₜ₋₁* ≥ 0}       0    0].
    ϕ1, ϕ2 = ϕs
    ϕ1_star, ϕ2_star = ϕ_stars
    Σ = []
    vals = [1, -1]
    for val in vals
        A = zeros(3, 3); A[2, 1] = 1; A[1, 2] = ϕ2_star; A[1, 3] = ϕ2
        A[1, 1] = ϕ1_star + ϕ1 * indicator(val, 0)
        A[3, 1] = indicator(val, 0)
        push!(Σ, A)
    end
    return Vector{Array{Float64, 2}}(Σ)
    # NOTE Σ has form s.t. Σ[1] is with yₜ₋₁ ≥ 0 and Σ[2] is with yₜ₋₁ < 0
end


function AR2_to_TAR(ϕs, ϕ_stars)
    # Returns the set of matrices 𝐀 from putting AR(2) model
    #        yₜ* = ϕ₁*yₜ₋₁* + ϕ₁yₜ₋₁ + ϕ₂*yₜ₋₂* + ϕ₂yₜ₋₂ +
    #        yₜ = max{yₜ*, 0}
    # in form
    #        𝐲ₜ* = 𝐀𝐲ₜ₋₁* + 𝛆ₜ
    # where   𝐲ₜ* = [yₜ*, yₜ₋₁*]; 𝛆ₜ = [ϵₜ, 0]
    # and     𝐀 = [ϕ₁*+ϕ₁𝐈{yₜ₋₁* ≥ 0}   ϕ₂*+ϕ₂𝐈{yₜ₋₂* ≥ 0}
    #                     1                     0        ].
    ϕ1, ϕ2 = ϕs
    ϕ1_star, ϕ2_star = ϕ_stars
    Σ = []
    vals = [1, -1]
    for val1 in vals
        for val2 in vals
            A = zeros(2, 2); A[2, 1] = 1
            A[1, 1] = ϕ1_star + ϕ1 * indicator(val1, 0)
            A[1, 2] = ϕ2_star + ϕ2 * indicator(val2, 0)
            push!(Σ, A)
        end
    end
    return Vector{Array{Float64, 2}}(Σ)
    # NOTE Σ is of form where Σ[1] is with yₜ₋₁,yₜ₋₂ ≥ 0, Σ[2] is with yₜ₋₁ ≥ 0 but yₜ₋₂ < 0, Σ[3] is with yₜ₋₁ < 0 but yₜ₋₂ ≥ 0 and Σ[4] is with yₜ₋₁,yₜ₋₂ < 0.
end

Σ = AR2_to_TAR([0.5,0.3], [0.2,0.1])
G = automaton_constructor(Σ)
s = discreteswitchedsystem(Σ, G, X)
sosbound_γ(s, 2)
using SwitchOnSafety, CSDP
s = discreteswitchedsystem(Σ, G)
soslyapb(s, 2, optimizer_constructor=CSDP.Optimizer)
