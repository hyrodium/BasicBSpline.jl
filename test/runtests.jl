using BasicBSpline
using IntervalSets
using Test

@testset "BasicBSpline.jl" begin
    # Write your own tests here.

    @testset "Knots" begin
        @test Knots(2,4,5) == Knots([2,4,5])
        @test zero(Knots) == Knots([])
        @test Knots(1:3) == Knots([3,2,1])
        @test Knots([-1,2,3]) + 2Knots([2,5]) == Knots([-1,2,2,2,3,5,5])
    end

    @testset "BSplineSpace" begin
        P1 = BSplineSpace(2,Knots([1,3,5,6,8,9]))
        @test BSplineSupport(2,P1) == 3..8
        @test dim(P1) == 3
        @test isproper(P1) == true

        P2 = BSplineSpace(1,Knots([1,3,3,3,8,9]))
        @test isproper(P2) == false
        @test properdim(P2) == 3

        i = 2
        k = Knots([5,12,13,13,14])
        p = 2
        P = 𝒫(p,k)
        @test BSplineSupport(P) == [5..13, 12..14]
        @test BSplineSupport(i,P) == 12..14
    end

end
