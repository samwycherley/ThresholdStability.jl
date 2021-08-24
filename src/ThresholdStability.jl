module ThresholdStability

using MathematicalSystems, SwitchOnSafety
using JuMP, SumOfSquares, DynamicPolynomials, MultivariatePolynomials
using CSDP
using Combinatorics, LinearAlgebra

import Reexport
Reexport.@reexport using HybridSystems

const AbstractSwitchedSystem = Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem, StateDepDiscreteSwitchedLinearSystem}

"""
    spectral_radius(A::AbstractMatrix)

Calculates the spectral radius of matrix `A`, using a utility from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
"""
spectral_radius(A::AbstractMatrix) = ρ(A)  # from SwitchOnSafety.jl

include("conversion.jl")

"""
    jsr(Σ::AbstractVector{<:AbstractMatrix})

Calculates an upper bound on the joint spectral radius of a set of matrices `Σ`, using a utility from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).

    jsr(s::AbstractSwitchedSystem)

Calculates an upper bound on the joint spectral radius for the switched linear system `s`, using a utility from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
"""
jsr(Σ::AbstractVector{<:AbstractMatrix}; d = 2, optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)) = soslyapb(discreteswitchedsystem(Σ), d, optimizer_constructor = optimizer_constructor)[2]
jsr(s::AbstractSwitchedSystem; d = 2, optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)) = soslyapb(unconstrained(unstatedep(s)), d, optimizer_constructor = optimizer_constructor)[2]

"""
    cjsr(Σ::AbstractVector{<:AbstractMatrix}, G::AbstractAutomaton)

Calculates an upper bound on the constrained joint spectral radius of a switching system `(Σ, G)`, using a utility from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
"""

cjsr(Σ::AbstractVector{<:AbstractMatrix}, G::AbstractAutomaton; d = 2, optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)) = soslyapb(discreteswitchedsystem(Σ, G), d, optimizer_constructor = optimizer_constructor)[2]
cjsr(s::AbstractSwitchedSystem; d = 2, optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)) = soslyapb(unstatedep(s), 2, optimizer_constructor = optimizer_constructor)[2]

export spectral_radius, cjsr, jsr

include("indicator.jl")
include("automaton_constructor.jl")

include("cksvar_functions.jl")

include("veronese_lift.jl")
include("midpoints.jl")
include("sdp_lyapunov_feasibility.jl")
include("sdp_gamma_search.jl")
include("sos_lyapunov_feasibility.jl")
include("sos_gamma_search.jl")
end
