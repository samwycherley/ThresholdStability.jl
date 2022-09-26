

import SwitchOnSafety: isinfeasible, isfeasible, isdecided, usestep
export sosbound_γ, sosbound_gamma


function sos_bisection_alg(s::StateDepDiscreteSwitchedLinearSystem, d, optimizer, lo, hi; 
        verbose=1, tol=1e-5, step=0.5)
    @assert hi > lo
    stuck_count = 0
    while hi - lo ≥ tol
        γ_mid = midpoint(hi, lo, step)
        status = soslyap_feasible_prog(γ_mid, s, d, optimizer, verbose=verbose-1)
        if verbose ≥ 1
            println("$γ_mid")
        end
        if !isdecided(status, true)
            if usestep(lo, hi)
                step *= 2
                continue
            end
            stuck_count += 1
            γ_midlo = min(midpoint(lo, γ_mid, step), shift(γ_mid, -tol/2))
            statuslo = soslyap_feasible_prog(γ_midlo, s, d, optimizer, verbose=verbose-1)
            if isdecided(statuslo, true)
                γ_mid = γ_midlo
                status = statuslo
            else
                γ_midhi = max(midpoint(γ_mid, hi, step), shift(γ_mid, tol/2))
                statushi = soslyap_feasible_prog(γ_midhi, s, d, optimizer, 
                    verbose=verbose-1)
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
        if stuck_count >= 5
            @warn("Making slow progress. Reporting upper bound with reduced accuracy.")
            return lo, hi
        end
    end
    return lo, hi
end


"""
    sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d)

Compute minimum value of ``γ`` for which the sum-of-squares program is feasible.

# Keywords
- `optimizer=...`: choose an SDP solver. The default solver is CSDP.
- `tol=...`: set tolerance. The default is `tol=1e-5`.
- `verbose=...`: set level of detail reported during calculation process. The default is 
`verbose=0`.
- `initstep=...`: set initial step size. The default value is `initstep=1.1`.
"""
function sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d; 
        optimizer=optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true), 
        tol=1e-5, verbose=0, initstep=1.1)
    lo_init, hi_init = 0, initstep; hi_stat = (MOI.INFEASIBLE, MOI.INFEASIBLE_POINT, 
        MOI.INFEASIBILITY_CERTIFICATE, "")
    while isinfeasible(hi_stat, true)
        hi_stat = soslyap_feasible_prog(hi_init, s, d, optimizer, verbose=verbose-1)
        if isinfeasible(hi_stat, true)
            lo_init = copy(hi_init)
            hi_init += initstep
        end
        initstep *= 2
    end
    lo, hi = sos_bisection_alg(s, d, optimizer, lo_init, hi_init, verbose=verbose)
    return hi
end

"""
    sosbound_gamma(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=nothing, tol=1e-5, 
    verbose=0, initstep=1.1)

Alias for [`sosbound_γ`](@ref).
"""
sosbound_gamma(s::StateDepDiscreteSwitchedLinearSystem, d; 
    optimizer=optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true), tol=1e-5, verbose=0, initstep=1.1) = sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d; 
        optimizer=optimizer, tol=tol, verbose=verbose, initstep=initstep)
