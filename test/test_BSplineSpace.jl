@testset "BSplineSpace" begin
    Random.seed!(42)

    P1 = BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9]))
    @test bsplinesupport(P1, 2) == 3..8
    @test dim(P1) == 3
    @test properdim(P1) == 3
    @test isproper(P1) == true

    P2 = BSplineSpace{1}(KnotVector([1, 3, 3, 3, 8, 9]))
    @test isproper(P2) == false
    @test properdim(P2) == 3

    i = 2
    k = KnotVector([5, 12, 13, 13, 14])
    p = 2
    P = BSplineSpace{p}(k)
    @test bsplinesupport(P) == [5..13, 12..14]
    @test bsplinesupport(P, i) == 12..14

    @test isproper(BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9])))
    @test !isproper(BSplineSpace{1}(KnotVector([1, 3, 3, 3, 8, 9])))

    @test dim(BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9]))) == 3

    P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
    P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
    P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
    @test P1 ⊆ P2
    @test P1 ⊆ P3
    @test P2 ⊈ P3

    P4 = BSplineSpace{1}(KnotVector([1, 2, 3, 4, 5]))
    P5 = BSplineSpace{2}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
    @test P4 ⊑ P4
    @test P4 ⊑ P5
    @test P5 ⊒ P4
    @test P5 ⋢ P4
    @test P4 ⋣ P5

    P4_ = BSplineSpace{degree(P4)-1}(knots(P4)[2:end-1])
    P5_ = BSplineSpace{degree(P5)-1}(knots(P5)[2:end-1])
    @test P4_ ⊑ P5_

end
