using DynamicPolynomials
using MathematicalSystems
import ThresholdStability: veroneselift

@testset "Veronese lifts" begin
    s = ContinuousIdentitySystem(4)
    @test veroneselift(s, 2) == s
    s = discreteswitchedsystem([[1. 2.; 1. 3.]], [[[1. 0.; 0. 1.], [1. -1.]]]).modes[1]
    @test veroneselift(s, 2).X[1] == I(3)
    s = discreteswitchedsystem([[1. 2.; 1. 3.]], [[[1. 0.; 0. 1.], [1. -1.]]]).modes[1]
    @test veroneselift(s, 2).X[2] ≈ [1 -√2 1]
    s = discreteswitchedsystem([[1. 2.; 1. 3.]], [[[1. 0.; 0. 1.], [1. -1.]]]).resetmaps[1]
    @test veroneselift(s, 2).A ≈ [1 √2*2 4; √2 5 √2*6; 1 √2*3 9]
    @polyvar x[1:2]
    @test veroneselift(x, 2) == [x[1]^2, √2x[1]*x[2], x[2]^2]
end
