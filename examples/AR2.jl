using ThresholdStability
using Plots, Distributions
using LaTeXStrings
pyplot()

E1, E2, E3, E4 = [1 0.; 0 1.], [1 0.; 0 -1.], [-1 0.; 0 1.], [-1 0.; 0 -1.]
D1 = zeros(1,2); D2, D3, D4 = copy(D1), copy(D1), copy(D1)
X = [[E1, D1], [E2, D2], [E3, D3], [E4, D4]]


function AR2_to_TAR(ϕs, ϕ_stars)
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

function AR2_to_companion(ϕs, ϕ_stars)
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

function simulate_AR2(y0, Σ, T, σ)  # for companion form
    y = zeros(3, T)
    y[:, 1] = y0
    for t in 1:T-1
        ϵ_t = [rand(Normal(0, σ)), 0., 0.]
        if y[1, t] ≥ 0.
            y[:, t+1] = Σ[1]*y[:, t] + ϵ_t
        else
            y[:, t+1] = Σ[2]*y[:, t] + ϵ_t
        end
    end
    return y
end

function plot_AR2(y0, Σ, T, σ; N = 20, row=1)
    # function to plot the latent variable in a censored/kinked AR(2) model
    ys = []
    ens_means = zeros(T)
    for i in 1:N
        y = simulate_AR2(y0, Σ, T, σ)
        y = y[row,:]
        push!(ys, y)
        ens_means .+= y
    end
    ens_means ./= N  # ensemble means

    E_y = simulate_AR2(y0, Σ, T, 0.)  # calculating deterministic results
    E_y = E_y[row,:]

    plot(ys, color = :grey, alpha = 0.1, label = "")
    plot!(ens_means, color = :grey, linewidth = 2, label = "Ensemble mean")
    plot!(E_y, color = :blue, linewidth = 2, linestyle = :dash, label = "Deterministic")
    plot!(xlabel=L"t", legend=:topright)
end

# ϕ₁ = 0.4; ϕ₁* = 0.2; ϕ₂ = 0.2; ϕ₂* = 0.1
Σ = AR2_to_TAR([0.4, 0.2], [0.2, 0.1])
G = automaton_constructor(Σ)
s = discreteswitchedsystem(Σ, G, X)
jsr(s)  # 0.925
cjsr(s)  # 0.925
sosbound_γ(s, 2)  # 0.925

Σ = AR2_to_companion([0.4, 0.2], [0.2, 0.1])
jsr(Σ)  # 0.925
plot_AR2(3*ones(3), Σ, 200, 1, row=1)
plot!(ylabel=L"y^*", yguidefontrotation=-90)


# ϕ₁ = 0.5; ϕ₁* = 0.5; ϕ₂ = -0.47; ϕ₂* = -0.5
Σ = AR2_to_TAR([0.5, -0.47], [0.5, -0.5])
G = automaton_constructor(Σ)
s = discreteswitchedsystem(Σ, G, X)
jsr(s)  # 1.105
cjsr(s)  # 1.001
sosbound_γ(s, 2)  # 0.985

Σ = AR2_to_companion([0.5, -0.47], [0.5, -0.5])
jsr(Σ)  # 1.003
plot_AR2(3*ones(3), Σ, 200, 1)
plot!(ylabel=L"y^*", yguidefontrotation=-90)

# ϕ₁ = 0.6; ϕ₁* = 0.6; ϕ₂ = -0.6; ϕ₂* = -0.6
Σ = AR2_to_TAR([0.6, -0.6], [0.6, -0.6])
G = automaton_constructor(Σ)
s = discreteswitchedsystem(Σ, G, X)
jsr(s)  # 1.245
cjsr(s)  # 1.118
sosbound_γ(s, 2)  # 1.095

Σ = AR2_to_companion([0.6, -0.6], [0.6, -0.6])
plot_AR2(3*ones(3), Σ, 200, 1)
plot!(ylabel=L"y^*", yguidefontrotation=-90)
