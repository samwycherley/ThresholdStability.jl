export CKSVAR_to_TAR, CKSVAR_to_companion, CKSVAR_to_companionFD

function split_by_var(C)
    # Splits π into πβ corresponding to πβ (all bar final row of π) and πβ corresponding to πbar*β (final row). I.e. π = [πβ; πβ].
    k = size(C,1)
    return C[1:k-1,:], C[k,:]'
end


function split_by_lag(C, nlags)
    # Given p lags, this splits π as per split_by_var and then separates πα΅’ into components πα΅’ΒΉ,β¦,πα΅’α΅ where πα΅’Κ² is the component of πα΅’ corresponding to the jth lag.
    if ndims(C) == 1
        k = size(C)
        n = 1
    else
        k, n = size(C)
    end
    @assert rem(n, nlags) == 0;
    m = Int(n / nlags)
    C1, C2 = split_by_var(C)
    C1s = []
    C2s = []
    for i in 1:nlags
        C1i = C1[:, (i-1)*m+1:i*m]
        C2i = C2[:, (i-1)*m+1:i*m]
        push!(C1s, C1i)
        push!(C2s, C2i)
    end
    return C1s, C2s
end


function atomize_C(C, nlags)
    # Subdivides components of π separated by lag into those corresponding to πββ and πββ.
    C1s, C2s = split_by_lag(C, nlags)
    C11s = []; C12s = []
    for (i, Ci) in enumerate(C1s)
        push!(C11s, C1s[i][:,1:end-1])
        push!(C12s, C1s[i][:,end])
    end
    C1s_atom = hcat(C11s, C12s)
    C21s = []; C22s = []
    for (i, Ci) in enumerate(C2s)
        push!(C21s, C2s[i][:,1:end-1])
        push!(C22s, C2s[i][:,end])
    end
    C2s_atom = hcat(C21s, C22s)
    return C1s_atom, C2s_atom
end

"""
    CKSVAR_to_TAR(C, Cstar, Ξ²tilde, nlags)

Returns pair `(Ξ£, X)` where `Ξ£` is a set of matrices corresponding to the CKSVAR model (with censored variable in levels) in TVAR form, and `X` is a set encoding state-space constraints.
"""
function CKSVAR_to_TAR(C, Cstar, Ξ²tilde, nlags)
    k,= size(C)
    C1s, C2s = atomize_C(C, nlags)
    Cs = reshape(hcat(C1s, C2s), (nlags, 2, 2))  # index notation is Cs[lag, rhsvariable, lhsvariable]
    Cstar1s, Cstar2s = split_by_lag(Cstar, nlags)
    Cstars = [Cstar1s Cstar2s]  # index notation is Cstars[lag, lhsvariable]

    vals = perms(nlags)

    Ξ£ = []
    Es = []; Ds = []; Xs = []
    for val in vals
        CTAR = zeros(k*nlags, k*nlags)
        CTAR[(k+1):end, 1:k*(nlags-1)] = I(k*(nlags-1))  # I
        CTAR[k+1:2k-1, k] = -Ξ²tilde * indicator(val[1], 0, indregion=:below)  # -Ξ²π{Yββββ < 0}

        CTAR2s = []
        CTAR2stars = []
        # lag = 1
        push!(CTAR2s, Cs[1, 1, 2])  # Cβββ
        push!(CTAR2stars, Cs[1, 2, 2] * indicator(val[1], 0) + Cstars[1, 2] - [Cs[1, 1, 2] * Ξ²tilde * indicator(val[1], 0, indregion=:below)])  # Cββ* + Cβββπ{Yβββα΅’ β₯ 0} - CβββΞ²π{Yβββα΅’ < 0}
        for lag in 2:nlags
            push!(CTAR2s, Cs[lag, 1, 2])  # Cββα΅’
            push!(CTAR2stars, Cstars[lag, 2] + Cs[lag, 2, 2] * indicator(val[lag], 0))  # Cβα΅’* + Cββα΅’π{Yβββα΅’ β₯ 0}
        end

        CTAR1s = []
        CTAR1stars = []
        # lag = 1
        push!(CTAR1s, Cs[1, 1, 1])  # Cβββ
        push!(CTAR1stars, Cstars[1, 1] + Cs[1, 2, 1] * indicator(val[1], 0) - Cs[1, 1, 1] * Ξ²tilde * indicator(val[1], 0 , indregion=:below))  # Cββ* + Cβββπ{Yβββα΅’* β₯ 0} - CβββΞ²π{Yβββα΅’* < 0}
        for lag in 2:nlags
            push!(CTAR1s, Cs[lag,1,1])  # Cββα΅’
            push!(CTAR1stars, Cstars[lag, 1] + Cs[lag, 2, 1] * indicator(val[lag], 0))  # Cβα΅’* + Cββα΅’π{Yβββα΅’* β₯ 0}
        end

        for i in 1:nlags
            CTAR[1:k, k*(i-1)+1:i*k] = hcat(vcat(CTAR1s[i], CTAR2s[i]), vcat(CTAR1stars[i], CTAR2stars[i]))
        end
        push!(Ξ£, CTAR)

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

    Vector{Array{Float64, 2}}(Ξ£), Vector(Xs)
    # NOTE the ordering of the matrices in Ξ£
end


"""
    CKSVAR_to_companion(C, Cstar, Ξ²tilde, nlags)

Converts CKSVAR model (with censored variable in levels) into companion form.
"""
function CKSVAR_to_companion(C, Cstar, Ξ²tilde, nlags)
    k = size(C,1)
    C1s, C2s = atomize_C(C, nlags)
    Cs = reshape(hcat(C1s, C2s), (nlags, 2, 2))  # index notation is Cs[lag, rhsvariable, lhsvariable]
    Cstar1s, Cstar2s = split_by_lag(Cstar, nlags)
    Cstars = [Cstar1s Cstar2s]  # index notation is Cstars[lag, lhsvariable]
    K = k + (nlags-1) * (k+1)
    Ξ£ = []
    for val in [1, -1]
        A = zeros(K, K)
        A[1:k-1, 1:k-1] = Cs[1, 1, 1]  # Cβββ
        A[k, 1:k-1] = Cs[1, 1, 2]  # Cβββ
        A[1:k-1, k] = Cstars[1, 1] + Cs[1, 2, 1]*indicator(val, 0, indregion=:above) - Cs[1, 1, 1]*Ξ²tilde*indicator(val, 0, indregion=:below)  # Cββ* + Cβββπ{Yββββ* β₯ 0} - CβββΞ²π{Yββββ* < 0}
        A[k,k] = Cstars[1, 2][1] + Cs[1, 2, 2][1]*indicator(val, 0, indregion=:above) - Cs[1, 1, 2]*Ξ²tilde*indicator(val, 0, indregion=:below)  # Cββ* + Cβββπ{Yββββ β₯ 0} - Cβββπ{Yββββ* < 0}
        for lag in 2:nlags
            C_lag = vcat(hcat(Cs[lag, 1, 1], Cs[lag, 2, 1], Cstars[lag, 1]), hcat(Cs[lag, 1, 2], Cs[lag, 2, 2], Cstars[lag, 2]))  # [Cββα΅’ Cββα΅’ Cβα΅’*; Cββα΅’ Cββα΅’ Cβα΅’*] (ith lag)
            J = (lag - 1)*(k + 1)
            A[1:k, J:J+k] = C_lag
        end
        A[k+1:2k-1, 1:k-1] = I(k-1)  # I
        A[k+1:2k-1, k] = -Ξ²tilde * indicator(val, 0, indregion=:below)  # -Ξ²π{Yββββ < 0}
        A[2k, k] = indicator(val, 0, indregion=:above)  # π{Yββββ β₯ 0}
        A[2k+1, k] = 1
        if nlags > 2
            A[2k+2:K, k+1:K-k-1] = I(K - 2k - 1)  # I
        end
        push!(Ξ£, A)
    end
    Vector{Array{Float64, 2}}(Ξ£)
end


"""
    CKSVAR_to_companionFD(F, Fstar, Ξ²tilde, nlags; diff=true)

Converts CKSVAR model estimated with the censored variable entering in first differences into companion form. The default setting `diff=true` returns the companion form with the censored variable entering in first differences. To retrieve the companion form with censored variable entering in levels, set `diff=false`.
"""
function CKSVAR_to_companionFD(F, Fstar, Ξ²tilde, nlags; diff = true)
    k = size(F, 1)
    F1s, F2s = atomize_C(F, nlags)
    Fs = reshape(hcat(F1s, F2s), (nlags, 2, 2))  # index notation is Fs[lag, rhsvariable, lhsvariable]
    Fstar1s, Fstar2s = split_by_lag(Fstar, nlags-1)
    Fstars = [Fstar1s Fstar2s]  # index notation is Fstars[lag, lhsvariable]
    if diff == false
        K = k + (nlags-1) * (k+1)
        Ξ£ = []
        for val in [1, -1]
            A = zeros(K, K)
            A[1:k-1, 1:k-1] = Fs[1, 1, 1]  # Cβββ
            A[k, 1:k-1] = Fs[1, 1, 2]  # Cβββ
            A[1:k-1, k] = Fstars[1, 1] + Ξ²tilde + (Fs[1, 2, 1] - Ξ²tilde)*indicator(val, 0, indregion=:above) - Fs[1, 1, 1]*Ξ²tilde*indicator(val, 0, indregion=:below)  # Fββ* + Ξ² + (Fβββ - Ξ²)π{Yββββ* β₯ 0} - CβββΞ²π{Yββββ* < 0}
            A[k,k] = 1. + Fstars[1, 2][1] + Fs[1, 2, 2][1]*indicator(val, 0, indregion=:above) - Fs[1, 1, 2]*Ξ²tilde*indicator(val, 0, indregion=:below)  # 1 + Fββ* + Fβββπ{Yββββ β₯ 0} - CβββΞ²π{Yββββ* < 0}
            for lag in 2:nlags-1  # i = 2,β¦,p-1
                F_lag = vcat(hcat(Fs[lag,1,1], Fs[lag,2,1] - Fs[lag-1,2,1], Fstars[lag, 1] - Fstars[lag-1,1]), hcat(Fs[lag, 1, 2], Fs[lag,2,2] - Fs[lag-1,2,2], Fstars[lag, 2] - Fstars[lag-1,2]))
                J = (lag - 1)*(k + 1)
                A[1:k, J:J+k] = F_lag
            end
            F_lag = vcat(hcat(Fs[nlags,1,1], -Fs[nlags-1,2,1], -Fstars[nlags-1,1]), hcat(Fs[nlags, 1, 2], -Fs[nlags-1,2,2], -Fstars[nlags-1,2]))  # lag = p
            J = (nlags - 1)*(k + 1)
            A[1:k, J:J+k] = F_lag
            A[k+1:2k-1, 1:k-1] = I + zeros(k-1, k-1)  # I
            A[k+1:2k-1, k] = -Ξ²tilde * indicator(val, 0, indregion=:below)  # -Ξ²π{Yββββ < 0}
            A[2k, k] = indicator(val, 0, indregion=:above)  # π{Yββββ β₯ 0}
            A[2k+1, k] = 1
            if nlags > 2
                A[2k+2:K, k+1:K-k-1] = I(K-2k-1)  # I
            end
            push!(Ξ£, A)
        end
    end
        if diff == true
            K = nlags * k
            Ξ£ = []
            for Ξ΄ in [1., 0.]
                A = zeros(K, K)
                A[1:k-1, 1:k-1] = Fs[1, 1, 1]  # Cβββ
                A[k, 1:k-1] = Fs[1, 1, 2]  # Cβββ
                A[1:k-1, k] = Fstars[1, 1] + Fs[1, 1, 1]*Ξ²tilde*(Ξ΄ - 1) + Fs[1, 2, 1]*Ξ΄  # Fββ* + CβββΞ²(Ξ΄ - 1) + FβββΞ΄
                A[k,k] = 1. + Fstars[1, 2][1] + Fs[1, 1, 2]*Ξ²tilde*(Ξ΄ - 1) + Fs[1, 2, 2][1]*Ξ΄  # Fββ* + CβββΞ²(Ξ΄ - 1) + FβββΞ΄
                for lag in 2:nlags-1  # i = 2,β¦,p-1
                    F_lag = vcat(hcat(Fs[lag,1,1], Fs[lag,2,1], Fstars[lag, 1]), hcat(Fs[lag, 1, 2], Fs[lag,2,2], Fstars[lag, 2]))
                    J = (lag - 1)*(k + 1)
                    A[1:k, J:J+k] = F_lag
                end
                F_lag = vcat(Fs[nlags,1,1], Fs[nlags, 1, 2])  # lag = p
                J = (nlags - 1)*(k + 1)
                A[1:k, J:J+k-2] = F_lag
                A[k+1:2k-1, 1:k-1] = I + zeros(k-1, k-1)  # I
                A[k+1:2k-1, k] = Ξ²tilde*(Ξ΄ - 1)  # Ξ²(Ξ΄ - 1)
                A[2k, k] = Ξ΄  # Ξ΄
                if nlags > 2
                    A[2k+1, k] = 1
                    A[2k+2:K, k+1:K-k-1] = I(K-2k-1)  # I
                end
                push!(Ξ£, A)
            end
        end
    Vector{Array{Float64, 2}}(Ξ£)
end
