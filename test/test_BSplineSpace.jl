@testset "BSplineSpace" begin
    Random.seed!(42)

    P1 = BSplineSpace(2, Knots([1, 3, 5, 6, 8, 9]))
    @test bsplinesupport(P1, 2) == 3..8
    @test dim(P1) == 3
    @test properdim(P1) == 3
    @test isproper(P1) == true

    P2 = BSplineSpace(1, Knots([1, 3, 3, 3, 8, 9]))
    @test isproper(P2) == false
    @test properdim(P2) == 3

    i = 2
    k = Knots([5, 12, 13, 13, 14])
    p = 2
    P = BSplineSpace(p, k)
    @test bsplinesupport(P) == [5..13, 12..14]
    @test bsplinesupport(P, i) == 12..14

    @test isproper(BSplineSpace(2, Knots([1, 3, 5, 6, 8, 9])))
    @test !isproper(BSplineSpace(1, Knots([1, 3, 3, 3, 8, 9])))

    @test dim(BSplineSpace(2, Knots([1, 3, 5, 6, 8, 9]))) == 3

    P1 = BSplineSpace(1, Knots([1, 3, 5, 8]))
    P2 = BSplineSpace(1, Knots([1, 3, 5, 6, 8, 9]))
    P3 = BSplineSpace(2, Knots([1, 1, 3, 3, 5, 5, 8, 8]))
    @test P1 ⊆ P2
    @test P1 ⊆ P3
    @test P2 ⊈ P3

    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    A12 = changebasis(P1, P2)
    A13 = changebasis(P1, P3)
    Δ12 = [bsplinebasis(P1, i, t) - sum(A12[i, j] * bsplinebasis(P2, j, t) for j in 1:n2) for i in 1:n1, t in 8 * rand(10)]
    Δ13 = [bsplinebasis(P1, i, t) - sum(A13[i, j] * bsplinebasis(P3, j, t) for j in 1:n3) for i in 1:n1, t in 8 * rand(10)]
    @test norm(Δ12) < 1e-14
    @test norm(Δ13) < 1e-14

    P4 = BSplineSpace(1, Knots([1, 2, 3, 4, 5]))
    P5 = BSplineSpace(2, Knots([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
    @test P4 ⊑ P4
    @test P4 ⊑ P5
    @test P5 ⊒ P4
    @test P5 ⋢ P4
    @test P4 ⋣ P5

    P4_ = BSplineSpace(degree(P4) - 1, knots(P4)[2:end-1])
    P5_ = BSplineSpace(degree(P5) - 1, knots(P5)[2:end-1])
    @test P4_ ⊑ P5_

    n4, n5 = dim(P4), dim(P5)
    A45 = changebasis(P4, P5)
    Δ45 = [bsplinebasis(P4, i, t) - sum(A45[i, j] * bsplinebasis(P5, j, t) for j in 1:n5) for i in 1:n1, t in 2 * rand(10) .+ 2]
    @test norm(Δ45) < 1e-14
end
