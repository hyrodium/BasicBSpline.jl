@testset "ChangeBasis" begin
    ε = 1e-14
    Random.seed!(42)

    P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
    P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
    P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
    @test P1 ⊆ P2
    @test P1 ⊆ P3
    @test P2 ⊈ P3

    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    A12 = changebasis(P1, P2)
    A13 = changebasis(P1, P3)
    @test size(A12) == (n1, n2)
    @test size(A13) == (n1, n3)
    Δ12 = [bsplinebasis(P1, i, t) - sum(A12[i, j] * bsplinebasis(P2, j, t) for j in 1:n2) for i in 1:n1, t in 8 * rand(10)]
    Δ13 = [bsplinebasis(P1, i, t) - sum(A13[i, j] * bsplinebasis(P3, j, t) for j in 1:n3) for i in 1:n1, t in 8 * rand(10)]
    @test norm(Δ12) < ε
    @test norm(Δ13) < ε

    p4 = 1
    p5 = 2
    P4 = BSplineSpace{p4}(KnotVector([1, 2, 3, 4, 5]))
    P5 = BSplineSpace{p5}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
    @test P4 ⊑ P4
    @test P4 ⊑ P5
    @test P5 ⊒ P4
    @test P5 ⋢ P4
    @test P4 ⋣ P5

    n4, n5 = dim(P4), dim(P5)
    A45 = changebasis(P4, P5)
    @test size(A45) == (n4, n5)
    Δ45 = [bsplinebasis(P4, i, t) - sum(A45[i, j] * bsplinebasis(P5, j, t) for j in 1:n5) for i in 1:n1, t in 2 * rand(10) .+ 2]
    @test norm(Δ45) < ε

    P6 = BSplineSpace{p4-1}(knots(P4)[2:end-1])
    P7 = BSplineSpace{p5-1}(knots(P5)[2:end-1])
    @test P6 ⊑ P7
    @test P6 ⊈ P7

    n6, n7 = dim(P6), dim(P7)
    A67 = changebasis(P6, P7)
    @test size(A67) == (n6, n7)
    Δ67 = [bsplinebasis(P6, i, t) - sum(A67[i, j] * bsplinebasis(P7, j, t) for j in 1:n7) for i in 1:n1, t in 2 * rand(10) .+ 2]
    @test norm(Δ67) < ε
end
