@testset "FastBSplineSpace" begin
    P1 = FastBSplineSpace(2, Knots([1, 3, 5, 6, 8, 9]))
    @test bsplinesupport(2, P1) == 3..8
    @test dim(P1) == 3
    @test properdim(P1) == 3
    @test isproper(P1) == true

    P2 = FastBSplineSpace(1, Knots([1, 3, 3, 3, 8, 9]))
    @test isproper(P2) == false
    @test properdim(P2) == 3

    i = 2
    k = Knots([5, 12, 13, 13, 14])
    p = 2
    P = FastBSplineSpace(p, k)
    @test bsplinesupport(P) == [5..13, 12..14]
    @test bsplinesupport(i, P) == 12..14

    @test isproper(FastBSplineSpace(2, Knots([1, 3, 5, 6, 8, 9])))
    @test !isproper(FastBSplineSpace(1, Knots([1, 3, 3, 3, 8, 9])))

    @test dim(FastBSplineSpace(2, Knots([1, 3, 5, 6, 8, 9]))) == 3

    P1 = FastBSplineSpace(1, Knots([1, 3, 5, 8]))
    P2 = FastBSplineSpace(1, Knots([1, 3, 5, 6, 8, 9]))
    P3 = FastBSplineSpace(2, Knots([1, 1, 3, 3, 5, 5, 8, 8]))
    @test P1 ⊆ P2
    @test P1 ⊆ P3
    @test P2 ⊈ P3
end
