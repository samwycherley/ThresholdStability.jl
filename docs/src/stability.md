# Stability

## Motivation
Let ``\{\mathscr{X}_i\}_{i=1}^m`` be a partition of state space ``\mathscr{X}\subseteq\mathbb{R}^p`` s.t. each ``\mathscr{X}_i`` takes the form of a convex polyhedron,
```math
\mathscr{X}_i=\{x\in\mathbb{R}^p\mid E_ix\geq_* 0, D_ix=0\},
```
where the inequality ``\geq_*`` may be strict or weak.

The tools in [ThresholdStability.jl](https://github.com/samwycherley/ThresholdStability.jl) are intended to determine whether a discrete-time model of form
```math
\begin{aligned}x_{t+1}&=\Phi(x_t)x_t\\
&:=\sum_{\sigma=1}^mA_\sigma\mathbf{1}\{x_t\in\mathscr{X}_i\}x_t\end{aligned}
```
is asymptotically stable, where ``\mathbf{1}\{x\in\mathscr{X}_i\}`` is an indicator function and ``A_i`` are given ``p\times p`` matrices.

Sufficient conditions for asymptotic stability, in descending order of conservativeness, are that
- the _joint spectral radius_ (JSR) of the set ``\Sigma=\{A_1,\dots,A_m\}`` is strictly less than 1;
- the _constrained joint spectral radius_ (CJSR) of the constrained switched linear system ``(\Sigma, G)``, where ``G`` is an automaton describing the possible transitions between states, is strictly less than 1;
- the '_state-constrained joint spectral radius_' (SCJSR) of the state-dependent switched linear system ``(\Sigma, G, \{\mathscr{X}_i\})`` is strictly less than 1.

In general, exact computation of the JSR, CJSR and SCJSR are NP-hard, but can be efficiently approximated via semidefinite programming or sum-of-squares programming (see e.g. [Parillo and Jadbabaie, 2008](https://arxiv.org/abs/0712.2887)).

## JSR and CJSR
We provide functions to compute tight upper bounds on the JSR and CJSR and to compute the spectral radius of individual matrices, each based on utilities from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
```@docs
spectral_radius
jsr
cjsr
```

The `jsr` function applied to a state-dependent switched system `(\Sigma, G, X)` or a constrained switched system `(\Sigma, G)` ignores the state-space constraints and automaton and reports the JSR of `\Sigma`. Likewise, the `cjsr` function applied to a state-dependent switched system ignores the state-space constraints. Applying the `cjsr` function to an unconstrained system is equivalent to applying `jsr`.

For information on how to construct switched systems `s`, please see the [HybridSystems documentation](https://blegat.github.io/HybridSystems.jl/stable/lib/methods/#Switched-Systems-1) or see the [examples/src folder](https://github.com/samwycherley/ThresholdStability.jl/tree/master/examples/src).

## SCJSR
This package contains two methods to compute upper bounds on the SCJSR, one solving a semidefinite program (SDP) directly and another solving a sum of squares (SOS) program.

```@docs
sosbound_γ
sosbound_gamma
sdpbound_γ
sdpbound_gamma
```

Since each ``\mathscr{X}_i`` is a convex polyhedron described by two matrices, `E_i` and `D_i`, each state space constraints should be supplied as the matrix pair `[E_i, D_i]`. For example, suppose we have a model
```math
begin{aligned}
y_t^*&=\phi_1^*y_{t-1}^*+\phi_1y_{t-1}+\phi_2^*y_{t-1}^*+\phi_2y_{t-1}+\epsilon_t,\\
y_t&=\max\{y_t^*,0\}.
\end{aligned}
```
This can be described either by a 3-dimensional state vector ``(y_t^*, y_{t-1}^*, y_{t-1})`` and a set of two ``3\times3`` matrices, or by a 2-dimensional state vector ``(y_t^*, y_{t-1}^*)`` and a set of four ``2\times2`` matrices:
- In the former case, the state space constraints are encoded as
```jldoctest
E_1, E_2, E_3, E_4 = [1 0 0.; 0 1 0.], [1 0 0.; 0 -1 0.], [-1 0 0.; 0 1 0.], [-1 0 0.; 0 -1 0.]
D_1, D_3 = [0 1 -1.], [0 1 -1.]; D_2, D_4 = [0 0 1.], [0 0 1.]
X = [[E_1, D_1], [E_2, D_2], [E_3, D_3], [E_4, D_4]]

# output
4-element Vector{Vector{Matrix{Float64}}}:
 [[1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 1.0 -1.0]]
 [[1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0 0.0 1.0]]
 [[-1.0 0.0 0.0; 0.0 1.0 0.0], [0.0 1.0 -1.0]]
 [[-1.0 0.0 0.0; 0.0 -1.0 0.0], [0.0 0.0 1.0]]
```
The matrices `E_i` describe the state space constraints on ``(y_t^*,y_{t-1}^*)`` and the matrices `D_i` enforce the constraint that ``x_{t-1}=\max\{x_{t-1}^*,0\}``.
- In the latter case, the state space constraints are encoded as
```jldoctest
E_1, E_2, E_3, E_4 = [1 0.; 0 1.], [1 0.; 0 -1.], [-1 0.; 0 1.], [-1 0.; 0 -1.]
D_1 = zeros(1,2); D_2, D_3, D_4 = copy(D_1), copy(D_1), copy(D_1)
X = [[E_1, D_1], [E_2, D_2], [E_3, D_3], [E_4, D_4]]

# output
4-element Vector{Vector{Matrix{Float64}}}:
 [[1.0 0.0; 0.0 1.0], [0.0 0.0]]
 [[1.0 0.0; 0.0 -1.0], [0.0 0.0]]
 [[-1.0 0.0; 0.0 1.0], [0.0 0.0]]
 [[-1.0 0.0; 0.0 -1.0], [0.0 0.0]]
```
Here, each `D_i` is a matrix of zeros since there are no censored variables in the state vector.

## SCJSR calculation details
### SDP program
The upper bound on the SCJSR found by [`sdpbound_γ`](@ref) involves finding the minimum value of the scalar ``\gamma`` s.t. the semidefinite program
```math
\begin{aligned}
Q_i-E_i^TU_iE_i+D_i^TZ_iD_i&\succ0\qquad\text{for }i=1,\dots,m,\\
\gamma^2Q_i-A_\sigma^TQ_jA_\sigma-E_i^TU_{ij}E_i+D_i^TZ_{ij}D_i&\succeq0\qquad\text{for }(i,j,\sigma)\in G,
\end{aligned}
```
is feasible, where ``U_i\geq0``, ``U_{ij}\geq0``, ``Z_i\succ0`` and ``Z_{ij}\succ0`` are otherwise arbitrary. Here, ``\geq`` denotes element-wise inequality, ``A\succ0`` indicates ``A`` is positive definite and ``A\succeq0`` indicates ``A`` is positive semidefinite.

### SOS program
The upper bound on the SCJSR found using [`sosbound_γ`](@ref) is generally less conservative than that found by [`sdpbound_γ`](@ref); [`sosbound_γ`](@ref) is thus recommended.

As in [Parillo and Jadbabaie (2008)](https://arxiv.org/abs/0712.2887), let ``x^{[d]}`` denote the ``d``-lift of vector ``x`` and let ``A^{[d]}`` be defined as the matrix such that ``(Ax)^{[d]}=A^{[d]}x^{[d]}``. Let ``y=x^{[d]}``. Now, [`sosbound_γ`](@ref) involves finding ``\gamma`` that minimizes
```math
\begin{aligned}
y^TQ_iy-(E_i^{[d]}y)^TU_iE_i^{[d]}y+(D_i^{[d]}y)^TZ_iD_i^{[d]}y-y^Ty&\enskip\text{ is SOS}\qquad\text{for }i=1,\dots,m,\\
\gamma^{2d}y^TQ_iy-(A_\sigma^{[d]}y)^TQ_jA_\sigma^{[d]}y\qquad\qquad&\\
-(E_i^{[d]}y)^TU_{ij}^{[d]}E_i^{[d]}y + (D_i^{[d]}y)^TZ_{ij}D_i^{[d]}y &\enskip\text{ is SOS}\qquad\text{for }(i,j,\sigma)\in G,
\end{aligned}
```
where ``U_i\geq0``, ``U_{ij}\geq0``, ``Z_i\succ0`` and ``Z_{ij}\succ0``.
