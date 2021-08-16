using CSDP, JuMP
optim = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

E1, E2, E3, E4 = [1 0.; 0 1.], [1 0.; 0 -1.], [-1 0.; 0 1.], [-1 0.; 0 -1.]
D1 = [0. 0.]; D2, D3, D4 = copy(D1), copy(D1), copy(D1)
X = [[E1, D1], [E2, D2], [E3, D3], [E4, D4]]

Σ1 = [[1. -1.; 2. 1.], [1. -1.; 1. 1.], [3. -2.; 1. 0.], [-2. 0.; 1. -1.]]
G1 = automaton_constructor(Σ1)
s1 = discreteswitchedsystem(Σ1, G1, X)
Σ2 = [[1.34 -1.34; 1. 0.], [1.34 -0.67; 1. 0.], [0.67 -1.34; 1. 0.], [0.67 -0.67; 1. 0.]]
G2 = automaton_constructor(Σ2)
s2 = discreteswitchedsystem(Σ2, G2, X)
Σ3 = [[1. -5.; 2. 0.], [2. -3.; 1. 0.], [0 -1.; 0. 1.], [2. 3.; -1. -4.]]
G3 = automaton_constructor(Σ3)
s3 = discreteswitchedsystem(Σ3, G3, X)

γsdp1 = sdpbound_γ(s1, optimizer=optim)
γsos1 = sosbound_γ(s1, 2, optimizer=optim)
γsdp2 = sdpbound_γ(s2, optimizer=optim)
γsos2 = sosbound_γ(s2, 2, optimizer=optim)
γsdp3 = sdpbound_γ(s3, optimizer=optim)
γsos3 = sosbound_γ(s3, 2, optimizer=optim)

@testset "SDP program" begin
    @test (γsdp1 ≥ 2 && γsdp1 ≤ 2.1)
    @test (γsdp2 ≥ 1.1 && γsdp2 ≤ 1.2)
    @test (γsdp3 ≥ 3.1 && γsdp3 ≤ 3.2)
end

@testset "SOS program" begin
    @test (γsos1 ≥ 2 && γsos1 ≤ 2.1)
    @test (γsos2 ≥ 1.1 && γsos2 ≤ 1.2)
    @test (γsos3 ≥ 3.1 && γsos2 ≤ 3.2)
end
