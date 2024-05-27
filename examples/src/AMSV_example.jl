using ThresholdStability
using CSV, HTTP, DataFrames, DataFramesMeta
using LinearAlgebra

extract_est(parname, df) = @subset(df, in([parname]).(:Parameter)).:Estimate[1]

## Aruoba, Mlikota, Schorfheide & Villalvazo (2022)
# Extract Aruoba data
AMSV = CSV.read(HTTP.get("https://raw.githubusercontent.com/samwycherley/ThresholdStability.jl/master/examples/src/estimates/estimates_AMSV22.csv
").body, DataFrame, copycols=true)


# Construct relevant objects
# A(1)
ψ_π = extract_est("psi_pi", AMSV)
ψ_z = extract_est("psi_z", AMSV)
α_S_Δ = extract_est("alpha_S_Delta", AMSV)
β_D = extract_est("beta_D", AMSV)
α_S = α_S_Δ + β_D
γ_D = -extract_est("minus_gamma_D", AMSV)

A = Array{Float64}(I(3))
A[1,2] = -γ_D
A[2,1] = -ψ_z
A[2,3] = 1.
A[3,1] = -ψ_π
A[3,2] = -β_D
A[3,3] = -α_S
A_21 = A[2:end, 1]

# Φ_ϵ(1)
Φ_ϵ1 = inv(A)
Φ_ϵ1_11 = Φ_ϵ1[1,1]  # partition...
Φ_ϵ1_12 = Φ_ϵ1[1,2:end]
Φ_ϵ1_21 = Φ_ϵ1[2:end,1]
Φ_ϵ1_22 = Φ_ϵ1[2:end,2:end]

# B_{⋅1}
B_1 = zeros(12)
for i in 1:12
    B_1[i] = extract_est(string("B_$i","_1"), AMSV)
end

# Φ(1)
Φ_1 = zeros(12, 3)
for i in 1:12
    for j in 2:3
        Φ_1[i,j] = extract_est(string("Phi_$i","_$j"), AMSV)
    end
end

Φ_12 = Φ_1[:,2:end]

for i in 1:12
    Φ_1[i,1] = B_1[i] - dot(Φ_12[i,:], A_21)
end

# Φ_12_Δ
Φ_12_Δ = zeros(2)
for i in 1:2
    Φ_12_Δ[i] = extract_est(string("Phi_12_$i","_Delta"), AMSV)
end

# μ
μ = zeros(3)
for i in 1:3
    μ[i] = extract_est("mu_$i", AMSV)
end

# D_ii, λs, etc
D_ii = zeros(3)
for i in 1:3
    D_ii[i] = extract_est(string("D_$i","_$i"), AMSV)
end

λ_ϕ = extract_est("lambda_phi", AMSV)
λ_0 = extract_est("lambda_0", AMSV)
λ_l = extract_est("lambda_l", AMSV)

ρ_ς = zeros(3)
for i in 1:3
    ρ_ς[i] = extract_est("rho_varsigma_$i", AMSV)
end

ς = zeros(3)
for i in 1:3
    ς[i] = extract_est("varsigma_$i", AMSV)
end

# Construct Φ(0)
Φ_0 = similar(Φ_1)
Φ_02 = Φ_0[:,2:end]

denom = B_1[1] - dot(Φ_12[1,:], A_21)
Φ_Δ_coeffs = zeros(15)
Φ_Δ_coeffs[1] = 1
for i in 2:12
    Φ_Δ_coeffs[i] = (B_1[i] - dot(Φ_12[i,:], A_21))/denom
end
Φ_Δ_coeffs[13] = (1 - dot(Φ_ϵ1_12, A_21))/denom
Φ_Δ_coeffs[14:15] = -Φ_ϵ1_22 * A_21/denom

for i in 1:12
    Φ_02[i,:] = Φ_12[i,:] + Φ_Δ_coeffs[i] * Φ_12_Δ
end
Φ_0[:,2:end] = Φ_02

for i in 1:12
    Φ_0[i,1] = B_1[i] - dot(Φ_02[i,:], A_21)
end

# Construct Φ_ϵ(0)
Φ_ϵ0_12 = Φ_ϵ1_12 + Φ_Δ_coeffs[13]*Φ_12_Δ
Φ_ϵ0_22 = Φ_ϵ1_22 + Φ_Δ_coeffs[14:end]*Φ_12_Δ'

# VAR(4)...
# ϕ_0⁺y_t⁺ + ϕ_0⁻y_t⁻ + Φ_0ˣx_t = sum_lags(...) + ϵ_t
# here we have y_1t*A_11 + y_2t'A_21 = x_t'B + ϵ_1t

k = 12
Ψ_0 = zeros(k, k+1) 
Ψ_0[1,1] = A[1,1]
Ψ_0[1,2] = A[1,1]
Ψ_0[2,3] = A_21[1]
Ψ_0[3, 4] = A_21[2]
Ψ_0[4:end,5:end] = I(k-3)

P, Q = render_canonical(Ψ_0)  # get P, Q [see (2.15) in DMW23]

function transformation_matrices(P, k)  # k is length of y_2t
    p = size(P, 1) - 1
    ψ_0yy⁺ = P[1,1]
    ψ_0yy⁻ = P[2,2]
    ψ_0xy⁺ = P[3:2+k,1]
    ψ_0xy⁻ = P[3:2+k,2]
    Ψ_0xx = P[3:2+k, 3:2+k]
    Ψ_0remainder = P[3+k:end, 3+k:end]

    Ψ_0⁺ = zeros(p,p)
    Ψ_0⁺[1,1] = ψ_0yy⁺
    Ψ_0⁺[2:1+k,1] = ψ_0xy⁺
    Ψ_0⁺[2:1+k,2:1+k] = Ψ_0xx
    Ψ_0⁺[2+k:end,2+k:end] = Ψ_0remainder

    Ψ_0⁻ = zeros(p,p)
    Ψ_0⁻[1,1] = ψ_0yy⁻
    Ψ_0⁻[2:1+k,1] = ψ_0xy⁻
    Ψ_0⁻[2:1+k,2:1+k] = Ψ_0xx
    Ψ_0⁻[2+k:end,2+k:end] = Ψ_0remainder    

    return Ψ_0⁺, Ψ_0⁻
end


# Turn Φ(1), Φ(0) into square matrices and convert to canonical
function canonical_transformation(Φ_0, Φ_1, P)
    p, q = size(Φ_0)
    Ψ_0⁺, Ψ_0⁻ = transformation_matrices(P, q-1)  # q is what we want
    Φ_0new = zeros(p,p)
    Φ_0new[1:q,1:end] = Φ_0'
    Φ_0new[q+1:end,1:end-q] = I(p-q)
    Φ_1new = zeros(p,p)
    Φ_1new[1:q,1:end] = Φ_1'
    Φ_1new[q+1:end,1:end-q] = I(p-q)
    Ψ_0⁻*Φ_0new, Ψ_0⁺*Φ_1new
end
Φ_0tilde, Φ_1tilde = canonical_transformation(Φ_0, Φ_1, P)
# These are the Φ(0), Φ(1) from AMSV with canonical transformation applied


import ThresholdStability: perms

# The representation here is with y_1t*s
function AMSV_to_TAR(Φ_1tilde, Φ_0tilde)
    nlags = 4
    k = 3
    vals = perms(nlags)

    Σ = []
    Es = []; Ds = []; Xs = []
    for val in vals
        if indicator(val[1], 0) == 1
            Φ_TAR = Φ_1tilde
        else
            Φ_TAR = Φ_0tilde
        end

        for lag in 2:nlags
            if indicator(val[lag], 0) == 0
                Φ_TAR[1:k,1+k*(lag-1)] = zeros(3)
            end
        end
        push!(Σ, Φ_TAR)

        E = zeros(nlags, k*nlags)
        D = zeros(1, k*nlags)
        for i in 1:nlags
            E[i, k*i] = val[i]
        end
        push!(Es, E)
        push!(Ds, D)
    end
    for i in 1:length(vals)
        X = Vector([Es[i], Ds[i]])
        push!(Xs, X)
    end

    Vector{Array{Float64, 2}}(Σ), Vector(Xs)
    # NOTE the ordering of the matrices in Σ
end

Σ, X = AMSV_to_TAR(Φ_1tilde, Φ_0tilde)
G = automaton_constructor(Σ)
s = discreteswitchedsystem(Σ, G, X)


AMSV_jsr = jsr(s)
AMSV_cjsr = cjsr(s)
AMSV_rjsr = rjsr(s)

using Tables
labels = ["jsr", "cjsr", "rjsr"]
vals = [AMSV_jsr, AMSV_cjsr, AMSV_rjsr]

results = Tables.table(hcat(labels, vals))

CSV.write("AMSV_results.csv", results, writeheader=false)  # save results