

import SwitchOnSafety: isinfeasible, isfeasible, isdecided, usestep
export sosbound_γ, sosbound_gamma

include("midpoints.jl")


function sos_bisection_alg(s::StateDepDiscreteSwitchedLinearSystem, d, optimizer, lo, hi; verbose=1, tol=1e-5, step=0.5)
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
                statushi = soslyap_feasible_prog(γ_midhi, s, d, optimizer, verbose=verbose-1)
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
    sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1)

Compute minimum value of ``\\gamma`` for which the sum-of-squares program is feasible.
"""

function sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1)
    if optimizer === nothing
        @warn("No optimizer supplied. Trying CSDP...")
        try
            optimizer = CSDP.Optimizer
        catch err
            @error("CSDP.Optimizer not found. Please supply an optimizer or install/rebuild CSDP.")
        end
    end
    lo_init, hi_init = 0, initstep; hi_stat = (MOI.INFEASIBLE, MOI.INFEASIBLE_POINT, MOI.INFEASIBILITY_CERTIFICATE, "")
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
    sosbound_gamma(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1)

Alias for [`sosbound_γ`](@ref).
"""


sosbound_gamma(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=nothing, tol=1e-5, verbose=0, initstep=1.1) = sosbound_γ(s::StateDepDiscreteSwitchedLinearSystem, d; optimizer=optimizer, tol=tol, verbose=verbose, initstep=initstep)
