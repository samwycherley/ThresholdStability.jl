# Replication code for numerical examples in Duffy, Mavroeidis and Wycherley,
# ``Stationarity with Occasionally Binding Constraints''

using ThresholdStability
using CSV, DataFrames, DataFramesMeta
using Parameters, Distributions

extract_est(parname, df) = @subset(df, in([parname]).(:Parameter)).:Estimate[1]

function render_canonical(Φ_0)
    ϕ_0yy_pos = Φ_0[1, 1]
    ϕ_0xy_pos = Φ_0[2:end, 1]
    ϕ_0yy_neg = Φ_0[1,2]
    ϕ_0xy_neg = Φ_0[2:end, 2]
    ϕ_0yx = Φ_0[1, 3:end]
    Φ_0xx = Φ_0[2:end, 3:end]
    d = size(Φ_0xx)
    Φ_0xx_inv = inv(Φ_0xx)

    ϕ_0yy_pos_canon = ϕ_0yy_pos - ϕ_0yx'*Φ_0xx_inv*ϕ_0xy_pos
    ϕ_0yy_neg_canon = ϕ_0yy_neg - ϕ_0yx'*Φ_0xx_inv*ϕ_0xy_neg

    P_inv = zeros(d[1]+2,d[2]+2)
    P_inv[1, 1] = ϕ_0yy_pos_canon
    P_inv[2, 2] = ϕ_0yy_neg_canon
    P_inv[3:end, 1] = ϕ_0xy_pos
    P_inv[3:end, 2] = ϕ_0xy_neg
    P_inv[3:end, 3:end] = Φ_0xx
    
    Q = Array{Float64}(I(d[2]+1))
    Q[1, 2:end] = -ϕ_0yx'*Φ_0xx_inv
    return (P = inv(P_inv), Q = Q)
end

## Example 2.1: monetary policy
MPModel = @with_kw (ψ = 0.9,                # natural rate persistence
                    γ = 1.5,                # Taylor rule coefficient
                    χ = 0.0,                # inflation persistence
                    θ = -0.5,               # responsiveness of π to r. N.b. θ<0
                    μ = 0.5,                # efficacy of unconventional MP
                    π_bar = 0.02,           # inflation target
                    r_bar = 0.04,           # mean natural rate
                    η_dist = Normal(0,0.1), # distribution for iid shocks η_t
                    ϵ_dist = Normal(0,0.1)  # distribution for iid shocks ϵ_t
)


function construct_mp_canonical(params)
    @unpack ψ, γ, μ, θ, χ = params
    κ_μ = 1/(1-μ*γ*θ)
    κ_1 = 1/(1-γ*θ)
    τ_μ = γ*θ*(1-μ)*κ_1
    Φ_pos = [ψ γ*(χ*κ_1-ψ)*κ_1; 0 χ*κ_1]
    Φ_neg = Φ_pos + [-χ*τ_μ*κ_μ 0; -χ*θ*(1-μ)*κ_μ 0]
    return [Φ_pos, Φ_neg]
end

# Toy parametrizations
# χ = 0.2, 0.5, 0.8
# μ = 0, 0.5, 1
# θ = -0.5, -1
# ψ = -0.5, 0, 0.5

χs = [0.2,0.7,0.99]
μs = [0.5,1]
θs = [-0.5,-1]
ψs = [0.1,0.5, 0.9]

outputs = []
for χ in χs
    for μ in μs
        for θ in θs
            for ψ in ψs
                model = MPModel(χ=χ,μ=μ,θ=θ,ψ=ψ)
                s = discreteswitchedsystem(construct_mp_canonical(model))
                ρ = jsr(s)
                push!(outputs, [(χ=χ,μ=μ,θ=θ,ψ=ψ), ρ])
            end
        end
    end
end
outputs[16]
# reported values correspond to indices 1,2,3, 16,17,18, 25,26,27

# Final case
model = MPModel(χ=0.9999, μ=0.5,ψ=0.999999,θ=-10)
@show jsr(discreteswitchedsystem(construct_mp_canonical(model)))
@show ρ.(construct_mp_canonical(model))


## Example 2.2: dynamic Tobit
E1, E2, E3, E4 = [1 0.; 0 1.], [1 0.; 0 -1.], [-1 0.; 0 1.], [-1 0.; 0 -1.]
D1 = zeros(1,2); D2, D3, D4 = copy(D1), copy(D1), copy(D1)

DynamicTobit = @with_kw (c = 0,                                         # constant
                         k = 2,                                         # number of lags
                         ϕ_pos = [0.9, 0.2],                            # vector ϕ_+
                         ϕ_neg = [0.5, 0.5],                            # vector ϕ_-
                         X = [[E1, D1], [E2, D2], [E3, D3], [E4, D4]]   # state space char
)


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
    # NOTE Σ is of form where Σ[1] is with yₜ₋₁,yₜ₋₂ ≥ 0, Σ[2] is with yₜ₋₁ ≥ 0 and yₜ₋₂ < 0, 
    # Σ[3] is with yₜ₋₁ < 0 and yₜ₋₂ ≥ 0 and Σ[4] is with yₜ₋₁,yₜ₋₂ < 0.
end


function construct_dynamic_tobit2(params)
    @unpack ϕ_pos, ϕ_neg, k = params
    @assert k == 2
    ϕ = ϕ_pos - ϕ_neg
    return AR2_to_TAR(ϕ, ϕ_neg)
end


function construct_dynamic_tobit2_system(params)
    @unpack X, k = params
    @assert k == 2
    Σ = construct_dynamic_tobit2(params)
    return discreteswitchedsystem(Σ, automaton_constructor(Σ), X)
end


# Parametrizations
# ϕ_pos = (.6, .3), ϕ_neg = (.2, .1)
dt_model = DynamicTobit(ϕ_pos = [.6, .3], ϕ_neg = [.2, .1])
s1 = construct_dynamic_tobit2_system(dt_model)
@show jsr(s1)
@show cjsr(s1)
@show rjsr(s1)

# ϕ_pos = (.6, .4), ϕ_neg = (.3, .1)
dt_model = DynamicTobit(ϕ_pos = [.6, .4], ϕ_neg = [.3, .1])
s2 = construct_dynamic_tobit2_system(dt_model)
@show jsr(s2)
@show cjsr(s2)
@show rjsr(s2)

# ϕ_pos = (.7, -.1), ϕ_neg = (.2, 0.)
dt_model = DynamicTobit(ϕ_pos = [.7, -.1], ϕ_neg = [.2, 0.])
s3 = construct_dynamic_tobit2_system(dt_model)
@show jsr(s3)
@show cjsr(s3)
@show rjsr(s3)

# ϕ_pos = (1.2, -1.2), ϕ_neg = (.6, -.6)
dt_model = DynamicTobit(ϕ_pos = [1.2, -1.2], ϕ_neg = [.6, -.6])
s4 = construct_dynamic_tobit2_system(dt_model)
@show jsr(s4)
@show cjsr(s4)
@show rjsr(s4)

# ϕ_pos = (1, -.97), ϕ_neg = (.5, -.5)
dt_model = DynamicTobit(ϕ_pos = [1, -.97], ϕ_neg = [.5, -.5])
s5 = construct_dynamic_tobit2_system(dt_model)
@show jsr(s5)
@show cjsr(s5)
@show rjsr(s5)

## Example 3.1
Φ_neg = [0 -1 0.36; 0 0.39 -1.33; 0 0.71 0.03]
Φ_pos = Φ_neg + [-1.37 0 0; 0.79 0 0; 0.76 0 0]

Σ = [Φ_neg, Φ_pos]
s = discreteswitchedsystem(Σ)
@show jsr(s)


## Ikeda, Li, Mavroeidis & Zanetti (2022): Japan estimates
JP = CSV.read("./estimates_JP.csv", DataFrame, copycols=true)  # Japan data from CSV

JPβtilde = zeros(2)
for i in 1:2
    JPβtilde[i] = extract_est("tildebeta_$i", JP)
end

JPCstar = zeros(3, 2)
for i in 1:3
    for j in 1:2
        JPCstar[i, j] = extract_est("Eq.$i lSR_$j", JP)
    end
end

Cbars = []
for j in 1:2
    Cbar_j = zeros(3, 3)
    for i in 1:3
        Cbar_j[i, 1] = extract_est("Eq.$i INFL_$j", JP)
        Cbar_j[i, 2] = extract_est("Eq.$i YGAP_BOJ_$j", JP)
        Cbar_j[i, 3] = extract_est("Eq.$i SR_$j", JP)
    end
    push!(Cbars, Cbar_j)
end
JPCbar = hcat(Cbars[1], Cbars[2])

JPC = copy(JPCbar); k = size(JPC, 1)  # convert C̄ to C
for i in 1:2
    JPC[:, i*k] -= JPCstar[:, i]
end

Σ_JP, X_JP = CKSVAR_to_TAR(JPC, JPCstar, JPβtilde, 2)

G_JP = automaton_constructor(Σ_JP)
s_JP = discreteswitchedsystem(Σ_JP, G_JP, X_JP)

@show jsr(s_JP)
@show cjsr(s_JP)
@show sosbound_γ(s_JP, 2)

