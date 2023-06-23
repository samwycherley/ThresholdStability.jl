using HybridSystems, MathematicalSystems

function unstatedep(s::AbstractSwitchedSystem)
    if typeof(s) <: StateDepDiscreteSwitchedLinearSystem
        n = nstates(s.automaton)
        oldmodes = s.modes
        modes = Vector{ContinuousIdentitySystem}(undef, n)
        for (i, mode) in enumerate(oldmodes)
            modes[i] = ContinuousIdentitySystem(mode.statedim)
        end
        G = s.automaton
        oldrmaps = s.resetmaps
        rmaps = Vector{LinearMap}(undef, n)
        for (i, map) in enumerate(oldrmaps)
            rmaps[i] = LinearMap(map.A)
        end
        sw = s.switchings
        s = HybridSystem(G, modes, rmaps, sw)
    end
    s
end


function unconstrained(s::AbstractSwitchedSystem)
    if typeof(s) <: ConstrainedDiscreteSwitchedLinearSystem
        s = HybridSystem(OneStateAutomaton(s.automaton.nt), s.modes, s.resetmaps, 
            s.switchings)
    elseif typeof(s) <: StateDepDiscreteSwitchedLinearSystem
        n = nstates(s.automaton)
        G = GraphAutomaton(n)
        for σ in 1:n
            for σnext in 1:n
                add_transition!(G, σ, σnext, σ)
            end
        end
        s = HybridSystem(G, s.modes, s.resetmaps, s.switchings)
    end
    s
end
