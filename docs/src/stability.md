# Stability

## Motivation
Let ``\{\mathscr{X}_\sigma\}_{\sigma=1}^m`` be a partition of state space ``\mathscr{X}\subseteq\mathbb{R}^p`` s.t. each ``\mathscr{X}_\sigma`` takes the form of a convex polyhedron,
```math
\mathscr{X}_\sigma=\{x\in\mathbb{R}^p\mid E_\sigma x\geq_* 0, D_\sigma x=0\},
```
where the inequality ``\geq_*`` may be strict or weak.

The tools in [ThresholdStability.jl](https://github.com/samwycherley/ThresholdStability.jl) are intended to determine whether a discrete-time model of form
```math
\begin{align}x_{t+1}&=\Phi(x_t)x_t\\
&:=\sum_{\sigma=1}^mA_\sigma\mathbf{1}\{x_t\in\mathscr{X}_\sigma\}x_t\end{align}
```
is asymptotically stable, where ``\mathbf{1}\{x\in\mathscr{X}_\sigma\}`` is an indicator function and ``A_\sigma`` are given ``p\times p`` matrices.

Sufficient conditions for asymptotic stability, in descending order of conservativeness, are that
- the _joint spectral radius_ (JSR) of the set ``\Sigma=\{A_1,\dots,A_m\}`` is strictly less than 1;
- the _constrained joint spectral radius_ (CJSR) of the constrained switched linear system ``(\Sigma, G)``, where ``G`` is an automaton describing the possible transitions between states, is strictly less than 1;
- the '_state-constrained joint spectral radius_' (SCJSR) of the state-dependent switched linear system ``(\Sigma, G, \{\mathscr{X}_\sigma\})`` is strictly less than 1.

In general, exact computation of the JSR, CJSR and SCJSR are NP-hard, but can be efficiently approximated via semidefinite programming or sum-of-squares programming (see e.g. [Parillo and Jadabaie, 2008](https://arxiv.org/abs/0712.2887)).

## JSR and CJSR
We provide functions to compute tight upper bounds on the JSR and CJSR and to compute the spectral radius of individual matrices, each based on utilities from [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl).
```@docs
spectral_radius
jsr
cjsr
```


The `jsr` function applied to a state-dependent switched system `(\Sigma, G, X)` or a constrained switched system `(\Sigma, G)` ignores the state-space constraints and automaton and reports the JSR of `\Sigma`. Likewise, the `cjsr` function applied to a state-dependent switched system ignores the state-space constraints.

## SCJSR

```@docs
sosbound_γ
sosbound_gamma
sdpbound_γ
sdpbound_gamma
```

## SCJSR calculation
