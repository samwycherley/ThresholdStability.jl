import ThresholdStability: unstatedep, unconstrained

Σ = [randn(2, 2) for _ in 1:4]
G = LightAutomaton(4); add_transition!(G,1,1,1); add_transition!(G,1,2,1); add_transition!(G,2,3,2); add_transition!(G, 3, 4, 3); add_transition!(G, 4, 1, 4)
X = [1,2,3,4]
s = discreteswitchedsystem(Σ, G, X)

@testset "Switched system type conversion" begin
    t = unstatedep(s)
    @test typeof(t) <: ConstrainedDiscreteSwitchedLinearSystem
    @test typeof(unconstrained(s)) <: StateDepDiscreteSwitchedLinearSystem
    @test typeof(unconstrained(t)) <: DiscreteSwitchedLinearSystem
end
