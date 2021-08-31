

import SwitchOnSafety: isinfeasible, isfeasible, isdecided, usestep
export sdpbound_γ, sdpbound_gamma


function sdpbisection_alg(s::StateDepDiscreteSwitchedLinearSystem, optimizer, lo, hi; verbose=1, tol=1e-5, step=0.5)
    @assert hi > lo
    while hi - lo ≥ tol
        γ_mid = midpoint(hi, lo, step)
        status = sdplyap_feasible_prog(γ_mid, s, optimizer, verbose=verbose-1)
        if verbose ≥ 1
            println("$γ_mid")
        end
        if !isdecided(status, true)
            if usestep(lo, hi)
                step *= 2
                continue
            end
            γ_midlo = min(midpoint(lo, γ_mid, step), shift(γ_mid, -tol/2))
            statuslo = sdplyap_feasible_prog(γ_midlo, s, optimizer, verbose=verbose-1)
            if isdecided(statuslo, true)
                γ_mid = γ_midlo
                status = statuslo
            else
                γ_midhi = max(midpoint(γ_mid, hi, step), shift(γ_mid, tol/2))
                statushi = sdplyap_feasible_prog(γ_midhi, s, optimizer, verbose=verbose-1)
                if isdecided(statushi, true)
                    γ_mid = γ_midhi
                    status = statushi
                end
            end
        end
        if isinfeasible(status, true)
            lo = γ_mid
        end
        if isfeasible(status, true)
            hi = γ_mid
        end
    end
    return lo, hi
end


"""
    sdpbound_γ(s::StateDepDiscreteSwitchedLinearSystem)

Compute minimum value of ``\\gamma`` for which the SDP is feasible.

For keywords, see [`sosbound_γ`](@ref)
"""
function sdpbound_γ(s::StateDepDiscreteSwitchedLinearSystem; optimizer=optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true), tol=1e-5, verbose=0, initstep=1.1)
    lo_init, hi_init = 0, initstep; hi_stat = (MOI.INFEASIBLE, MOI.INFEASIBLE_POINT, MOI.INFEASIBILITY_CERTIFICATE, "")
    while isinfeasible(hi_stat, true)
        hi_stat = sdplyap_feasible_prog(hi_init, s, optimizer, verbose=verbose-1)
        if isinfeasible(hi_stat, true)
            lo_init = copy(hi_init)
            hi_init += initstep
        end
        initstep *= 2
    end
    lo, hi = sdpbisection_alg(s, optimizer, lo_init, hi_init, verbose=verbose)
    return hi
end


"""
    sdpbound_gamma(s::StateDepDiscreteSwitchedLinearSystem; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1)

Alias for [`sdpbound_γ`](@ref).
"""
sdpbound_gamma(s::StateDepDiscreteSwitchedLinearSystem; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1) = bound_γ(s::StateDepDiscreteSwitchedLinearSystem; optimizer=optimizer, tol=tol, verbose=verbose, initstep=initstep)
