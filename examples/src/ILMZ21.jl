# # Testing the stability of a CKSVAR model

# Here, we test the stability of the CKSVAR model estimated for Japan in:
# - Ikeda, Li, Mavroeidis and Zanetti (2021), "[Testing the effectiveness of monetary policy in Japan and the United States](https://arxiv.org/abs/2012.15158)." Working paper.
# This paper is hereafter referred to as ILMZ21.
##
# In ILMZ21, the effectiveness of unconventional monetary policy in Japan is investigated via estimation of a CKSVAR(2) model in inflation, output gap and the Call Rate for Japan, over the period 1985q3-2019q1, with shadow call rate as the latent variable.
#
# Converting this model into TAR form, we determine upper bounds on the state-constrained joint spectral radius (SCJSR)
##
using DataFrames, DataFramesMeta, CSV, HTTP
using ThresholdStability

extract_est(parname, df) = @subset(df, in([parname]).(:Parameter)).:Estimate[1]

JP = CSV.read(HTTP.get("https://raw.githubusercontent.com/samwycherley/ThresholdStability.jl/master/examples/src/estimates/estimates_JP.csv").body, DataFrame, copycols=true)  # Japan data from CSV
##
##
# This loads the Japan estimates as a dataframe and defines a function to extract estimates. Next, we build the matrices:
##
##
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
##
# Write this as a state-dependent switched linear system:
##
##
G_JP = automaton_constructor(Σ_JP)
s_JP = discreteswitchedsystem(Σ_JP, G_JP, X_JP)
##
@show sosbound_γ(s_JP, 2)
##
# We see that the model estimated for Japan is stable