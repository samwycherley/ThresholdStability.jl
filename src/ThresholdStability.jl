module ThresholdStability

using MathematicalSystems, SwitchOnSafety
using JuMP, SumOfSquares, DynamicPolynomials, MultivariatePolynomials
using CSDP
using Combinatorics, LinearAlgebra

"""
    spectral_radius(A::AbstractMatrix)

Calculates the spectral radius of matrix `A`, using a utility from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
"""
spectral_radius(A::AbstractMatrix) = ρ(A)  # from SwitchOnSafety.jl
export spectral_radius

import Reexport
Reexport.@reexport using HybridSystems
Reexport.@reexport import SwitchOnSafety: soslyapb

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
