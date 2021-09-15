# # Testing the stability of CKSVAR models

# Here, we test the stability of the models estimated in:
# - Ikeda, Li, Mavroeidis and Zanetti (2021), "[Testing the effectiveness of monetary policy in Japan and the United States](https://arxiv.org/abs/2012.15158)." Working paper.
# This paper is hereafter referred to as [ILMZ21].
##
# # Japan
# In [ILMZ21], the effectiveness of unconventional monetary policy in Japan is investigated via estimation of a CKSVAR(2) model in inflation, output gap and the Call Rate for Japan, over the period 1985q3-2019q1, with shadow call rate as the latent variable.
#
# Converting this model into TAR form, we determine upper bounds on the state-constrained joint spectral radius (SCJSR)
##
using DataFrames, DataFramesMeta, CSV, HTTP
using ThresholdStability

extract_est(parname, df) = @subset(df, in([parname]).(:Parameter)).:Estimate[1]

JP = CSV.read(HTTP.get("https://raw.githubusercontent.com/samwycherley/ThresholdStability.jl/master/examples/src/estimates/estimates_JP.csv").body, DataFrame, copycols=true)  # Japan data from CSV

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
##
@show sosbound_γ(s_JP, 2)
##
# # United States
# For the United States, [ILMZ21] estimate a CKSVAR(3) model in inflation, output gap, and the Federal Funds Rate in first differences. The shadow Federal Funds Rate is the latent variable, and also enters in first differences.
#
# Because the Federal Funds Rate enters in first differences, the SCJSR is no longer appropriate, since `Σ` is now a continuum of matrices. Instead we rely on the joint spectral radius.
##
USfd = CSV.read(HTTP.get("https://raw.githubusercontent.com/samwycherley/ThresholdStability.jl/master/examples/src/estimates/estimates_US-fd.csv").body, DataFrame, copycols=true)  # US data from CSV

USfdβtilde = zeros(2)
for i in 1:2
    USfdβtilde[i] = extract_est("tildebeta_$i", USfd)
end

Cstar = zeros(3, 3)
for i in 1:3
    for j in 1:3
        Cstar[i, j] = extract_est("Eq.$i lFEDFUNDS_$j", USfd)
    end
end
USfdFstar = zeros(3, 2)
for i in 1:2
    for j in 1:i
        USfdFstar[1:2, i] += Cstar[1:2, j]
        USfdFstar[3, i] += Cstar[3, j]
    end
    USfdFstar[1:2, i] -= USfdβtilde
    USfdFstar[3, i] -= 1.
end

Cbars = []
for j in 1:3
    Cbar_j = zeros(3, 3)
    for i in 1:3
        Cbar_j[i, 1] = extract_est("Eq.$i INFL_$j", USfd)
        Cbar_j[i, 2] = extract_est("Eq.$i YGAP_$j", USfd)
        Cbar_j[i, 3] = extract_est("Eq.$i FEDFUNDS_$j", USfd)
    end
    push!(Cbars, Cbar_j)
end
Cs = []
for j in 1:3
    C_j = copy(Cbars[j])
    C_j[:, 3] -= Cstar[:, j]
    push!(Cs, C_j)
end
Fs = []
for i in 1:3
    F_i = zeros(3, 3)
    F_i[1:3, 1:2] = Cs[i][1:3, 1:2]
    for j in 1:i
        F_i[1:2, 3] += Cs[j][1:2, 3]
    end
    F_i[1:2, 3] += USfdβtilde
    for j in 1:i
        F_i[3, 3] += Cs[j][3, 3]
    end
    push!(Fs, F_i)
end
USfdF = hcat(Fs[1], Fs[2], Fs[3])

Σ_US = CKSVAR_to_companionFD(USfdF, USfdFstar, USfdβtilde, 3, diff=true)
##
@show jsr(Σ_US)
##
# Both models are stable.
