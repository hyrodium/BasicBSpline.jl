using BasicBSpline
using IntervalSets
using Test

@testset "BasicBSpline.jl" begin
    # Write your own tests here.

    @testset "Knots" begin
        @test zero(Knots)==Knots([])
        @test Knots(1:3)==Knots([3,2,1])
        @test Knots([-1,2,3])+2Knots([2,5])==Knots([-1,2,2,2,3,5,5])
        @test Knots([-1,2,3])+2Knots([2,5])==Knots([-1,2,2,2,3,5,5])
    end

    @testset "BSplineSpace" begin
        P = BSplineSpace(2,Knots([1,3,5,6,8,9]))
        @test BSplineSupport(2,P)==3..8
        @test dim(P)==3
    end

end
