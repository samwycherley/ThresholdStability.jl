

function constrain_sdp(model, s::ConstrainedContinuousIdentitySystem, Qi, Ui, Zi, nQ, nZ)
    Ei = s.X[1]; Di = s.X[2];
    @constraint(model, Qi - Ei'*Ui*Ei - Di'*Zi*Di - Array(I(nQ)) in PSDCone())
    @constraint(model, Zi - Array(I(nZ)) in PSDCone())
end


function constrain_sdp(model, s::ConstrainedLinearMap, Qi, Qj, Uij, Zij, nZ, γ)
    Ei = s.X[1]; Di = s.X[2]; Ai = s.A
    @constraint(model, γ^2*Qi - Ai'*Qj*Ai - Ei'*Uij*Ei - Di'*Zij*Di in PSDCone())
    @constraint(model, Zij - Array(I(nZ)) in PSDCone())
end


function sdplyap_feasible_prog(γ, s::StateDepDiscreteSwitchedLinearSystem, optimizer; verbose=0)
    modes=1:nstates(s)
    model = Model(optimizer)
    nA = size(s.resetmaps[1].A, 1)
    nE = size(s.resetmaps[1].X[1], 1)
    nD = size(s.resetmaps[1].X[2], 1)
    if verbose ≤ 0
        set_silent(model)
    end
    Q_vrefs = Dict(state => @variable(model, [i=1:nA, j=1:nA]) for state in modes)
    Ui_vrefs = Dict(state => @variable(model, [i=1:nE, j=1:nE], lower_bound=0) for state in modes)
    Zi_vrefs = Dict(state => @variable(model, [i=1:nD, j=1:nD]) for state in modes)

    for i in modes
        Qi = Q_vrefs[i]
        constrain_sdp(model, HybridSystems.mode(s, i), Qi, Ui_vrefs[i], Zi_vrefs[i], nA, nD)
        for t in out_transitions(s, i)
            j = target(s, t)
            Qj = Q_vrefs[j]
            Uij_ref = @variable(model, [i=1:nE, j=1:nE], lower_bound=0)
            Zij_ref = @variable(model, [i=1:nD, j=1:nD])
            constrain_sdp(model, resetmap(s, t),  Qi, Qj, Uij_ref, Zij_ref, nD, γ)
        end
    end


    if verbose ≥ 2
        print(model)
    end

    optimize!(model)

    if verbose ≥ 3
        print(model)
    end
    if verbose ≥ 2
        try
            @show JuMP.solution_summary(model, verbose=true)
        catch err
            println("solution_summary requires JuMP version 0.21.7 or later")
        end
    elseif verbose ≥ 1
        try
            @show JuMP.solution_summary(model)
        catch err
            println("solution_summary requires JuMP version 0.21.7 or later")
        end
    end

    status = (JuMP.termination_status(model), JuMP.primal_status(model), JuMP.dual_status(model), JuMP.raw_status(model))
end
