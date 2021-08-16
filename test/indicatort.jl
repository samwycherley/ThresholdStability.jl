
A = -5:0.2:5
B = -5:5
@testset "Indicator" begin
    @testset "above, nonstrict $y, $b" for y in A, b in B
        if y ≥ b
            @test indicator(y, b) == 1
        else
            @test indicator(y, b) == 0
        end
    end
    @testset "above, strict $y, $b" for y in A, b in B
        if y > b
            @test indicator(y, b, ineq=:strict) == 1
        else
            @test indicator(y, b, ineq=:strict) == 0
        end
    end
    @testset "below, nonstrict $y, $b" for y in A, b in B
        if y ≤ b
            @test indicator(y, b, indregion=:below) == 1
        else
            @test indicator(y, b, indregion=:below) == 0
        end
    end
    @testset "below, strict $y, $b" for y in A, b in B
        if y < b
            @test indicator(y, b, indregion=:below, ineq=:strict) == 1
        else
            @test indicator(y, b, indregion=:below, ineq=:strict) == 0
        end
    end
end
