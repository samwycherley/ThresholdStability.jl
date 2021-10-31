<meta name="google-site-verification" content="h3Hd4_8J-eewzmXcrLDC0Sa9Vp77aKmg6IabUXc9ObA" />
# ThresholdStability

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://samwycherley.github.io/ThresholdStability.jl/dev)
[![Build Status](https://github.com/samwycherley/ThresholdStability.jl/workflows/CI/badge.svg)](https://github.com/samwycherley/ThresholdStability.jl/actions)
[![Coverage](https://codecov.io/gh/samwycherley/ThresholdStability.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/samwycherley/ThresholdStability.jl)

A Julia package implementing techniques to determine the stability of discrete-time threshold vector autoregressive (TVAR) models such as censored and kinked vector autoregressive models (CKSVAR).

## Resources
- Documentation: <https://samwycherley.github.io/ThresholdStability.jl/dev/>.
- Examples: examples, including notebooks, are contained in the [examples/src folder](https://github.com/samwycherley/ThresholdStability.jl/tree/master/examples/src).
- R wrapper for this package: <https://github.com/samwycherley/thresholdr>
## Installation
ThresholdStability.jl requires Julia v1.0 or later. Julia can be downloaded [here](https://julialang.org/downloads/).

To install, use
```julia
Pkg.add("ThresholdStability")
```
or from the REPL, type
```julia
] add ThresholdStability
```

To test the package once installed, use
```julia
Pkg.test("ThresholdStability")
```
or from the REPL,
```julia
] test ThresholdStability
```
Testing after installation is recommended.

## Dependencies
The package relies on a number of other packages, most importantly
- [HybridSystems.jl](https://github.com/blegat/HybridSystems.jl)
- [SwitchOnSafety.jl](https://github.com/blegat/SwitchOnSafety.jl)
- [MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl)
- [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl)

All dependent packages are automatically installed when installing ThresholdStability.

For troubleshooting, consulting the documentation of these packages may be helpful, in addition to the [documentation for ThresholdStability.jl](https://samwycherley.github.io/ThresholdStability.jl/dev/).
