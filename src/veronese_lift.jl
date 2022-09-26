import SwitchOnSafety: veroneselift
import Base: copy


copy(Lmap::LinearMap) = LinearMap(copy(Lmap.A))

copy(Lmap::ConstrainedLinearMap) = ConstrainedLinearMap(copy(Lmap.A), copy(Lmap.X))

copy(s::ConstrainedContinuousIdentitySystem) = 
    ConstrainedContinuousIdentitySystem(s.statedim, copy(s.X))

copy(s::AbstractSwitchedSystem) = HybridSystem(s.automaton, copy.(s.modes), 
    copy.(s.resetmaps), s.switchings)


veroneselift(Lmap::LinearMap, d) = LinearMap(veroneselift(Lmap.A, d))

function veroneselift(Lmap::ConstrainedLinearMap, d)
    A = veroneselift(Lmap.A, d)
    X = Lmap.X
    if typeof(X) <: AbstractArray
        for (i, E) in enumerate(X)
            X[i] = veroneselift(E, d)
        end
    end
    ConstrainedLinearMap(A, X)
end


veroneselift(s::ContinuousIdentitySystem, d) = s

function veroneselift(s::ConstrainedContinuousIdentitySystem, d)
    X = copy(s.X)
    if typeof(X) <: AbstractArray
        for (i, E) in enumerate(s.X)
            X[i] = veroneselift(E, d)
        end
    end
    ConstrainedContinuousIdentitySystem(s.statedim, X)
end


function veroneselift!(s::AbstractSwitchedSystem, d)
    for (i, mode) in enumerate(s.modes)
        s.modes[i] = veroneselift(mode, d)
    end
    for (i, Rmap) in enumerate(s.resetmaps)
        s.resetmaps[i] = veroneselift(Rmap, d)
    end
end


function veroneselift(s::StateDepDiscreteSwitchedLinearSystem, d)
    t = copy(s)
    veroneselift!(t, d)
    t
end

const AbstractPolyVec = AbstractVector{<:MultivariatePolynomials.AbstractVariable}

function veroneselift(x::AbstractPolyVec, d)
    n = length(x)
    N = binomial(n + d - 1, d)
    X = monomials(x, d)
    df = factorial(d)
    scaling(n) = sqrt(df / prod(factorial, exponents(n)))
    scales = Float64[scaling(n) for n in X]
    return Vector(terms(scales ⋅ X))
end
