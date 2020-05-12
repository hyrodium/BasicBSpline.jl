using BasicBSpline
using IntervalSets
using Test

@testset "BasicBSpline.jl" begin
    # Write your own tests here.

    @testset "Knots" begin
        @test Knots([1,2,3]) == Knots(1:3) == Knots(1,2,3)
        @test Knots(2,4,5) == Knots([2,4,5])
        @test zero(Knots) == Knots([])
        @test Knots(1:3) == Knots([3,2,1])

        k = Knots([4,5,6])
        @test ♯(k) == length(k) == 3

        @test Knots([-1,2,3]) + 2Knots([2,5]) == Knots([-1,2,2,2,3,5,5])
        Knots([1,2,3]) + Knots([2,4,5]) == Knots([1,2,2,3,5])
        2 * Knots([2,3]) == Knots([2,2,3,3])

        unique(Knots([1,2,2,3])) == Knots([1,2,3])

        k = Knots([1,2,2,3])
        @test 𝔫(k,0.3) == 0
        @test 𝔫(k,1.0) == 1
        @test 𝔫(k,2.0) == 2
    end

    @testset "BSplineSpace" begin
        P1 = BSplineSpace(2,Knots([1,3,5,6,8,9]))
        @test bsplinesupport(2,P1) == 3..8
        @test dim(P1) == 3
        @test properdim(P1) == 3
        @test isproper(P1) == true

        P2 = BSplineSpace(1,Knots([1,3,3,3,8,9]))
        @test isproper(P2) == false
        @test properdim(P2) == 3

        i = 2
        k = Knots([5,12,13,13,14])
        p = 2
        P = 𝒫(p,k)
        @test bsplinesupport(P) == [5..13, 12..14]
        @test bsplinesupport(i,P) == 12..14

        @test isproper(𝒫(2,Knots([1,3,5,6,8,9])))
        @test !isproper(𝒫(1,Knots([1,3,3,3,8,9])))

        @test dim(𝒫(2,Knots([1,3,5,6,8,9]))) == 3

        P1 = 𝒫(1,Knots([1,3,5,8]))
        P2 = 𝒫(1,Knots([1,3,5,6,8,9]))
        P3 = 𝒫(2,Knots([1,1,3,3,5,5,8,8]))
        @test P1 ⊆ P2
        @test P1 ⊆ P3
        @test P2 ⊈ P3
    end

    @testset "FastBSplineSpace" begin
        P1 = FastBSplineSpace(2,Knots([1,3,5,6,8,9]))
        @test bsplinesupport(2,P1) == 3..8
        @test dim(P1) == 3
        @test properdim(P1) == 3
        @test isproper(P1) == true

        P2 = FastBSplineSpace(1,Knots([1,3,3,3,8,9]))
        @test isproper(P2) == false
        @test properdim(P2) == 3

        # i = 2
        # k = Knots([5,12,13,13,14])
        # p = 2
        # P = f𝒫(p,k)
        # @test bsplinesupport(P) == [5..13, 12..14]
        # @test bsplinesupport(i,P) == 12..14
        #
        # @test isproper(𝒫(2,Knots([1,3,5,6,8,9])))
        # @test !isproper(𝒫(1,Knots([1,3,3,3,8,9])))
        #
        # @test dim(𝒫(2,Knots([1,3,5,6,8,9]))) == 3
        #
        # P1 = 𝒫(1,Knots([1,3,5,8]))
        # P2 = 𝒫(1,Knots([1,3,5,6,8,9]))
        # P3 = 𝒫(2,Knots([1,1,3,3,5,5,8,8]))
        # @test P1 ⊆ P2
        # @test P1 ⊆ P3
        # @test P2 ⊈ P3
    end

    @testset "refinement" begin
        P1 = 𝒫(1,Knots([0,0,1,1]))
        P2 = 𝒫(1,Knots([1,1,2,3,3]))
        n1 = dim(P1) # 2
        n2 = dim(P2) # 3
        𝒂 = [[i, j] for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂ array.
        M = BSplineManifold([P1, P2], 𝒂)

        P1′ = 𝒫(2,Knots([0,0,0,1,1,1]))
        P2′ = 𝒫(1,Knots([1,1,2,1.45,3,3]))

        @test P1 ⊆ P1′
        @test P2 ⊆ P2′

        M′ = refinement(M, [P1′, P2′])
        t = [0.82,1.8]
        @test mapping(M, t) ≈ mapping(M′, t)
    end
end
