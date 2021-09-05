import ThresholdStability: perms

@testset "Ordering" begin
    @test perms(1) == [[1.],[-1.]]
    @test perms(2) == [[1.,1.],[1.,-1.],[-1.,1.],[-1.,-1.]]
    @test perms(3) == [[1.,1.,1.],[1.,1.,-1.],[1.,-1.,1.],[-1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.],[-1.,-1.,-1.]]
end

function Base.:≈(X::T, Y::T) where T <: LightAutomaton
    if X.G == Y.G && X.Σ == Y.Σ && X.nt == Y.nt && X.next_id == Y.next_id
        return true
    end
    false
end


Σ = [randn(2,2) for _ in 1:4]
G = LightAutomaton(4)
add_transition!(G, 1, 1, 1)
add_transition!(G, 1, 3, 1)
add_transition!(G, 2, 1, 2)
add_transition!(G, 2, 3, 2)
add_transition!(G, 3, 2, 3)
add_transition!(G, 3, 4, 3)
add_transition!(G, 4, 2, 4)
add_transition!(G, 4, 4, 4)


@testset "Automata" begin
    @test G ≈ automaton_constructor(Σ)
    for n in 2:7
        Σn = [randn(2,2) for _ in 1:2^n]
        @test automaton_constructor(Σn) ≈ automaton_constructor(Σn, perms(n))
    end
    @test islight(G) == true
    @test islight(OneStateAutomaton(2)) == false
end
