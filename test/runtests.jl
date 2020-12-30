using BasicBSpline
using IntervalSets
using LinearAlgebra
using Test
using Random

ε = 1.0e-8

@testset "BasicBSpline.jl" begin

    @testset "Knots" begin
        @test Knots([1, 2, 3]) == Knots(1:3) == Knots(1, 2, 3)
        @test Knots(2, 4, 5) == Knots([2, 4, 5])
        @test zero(Knots) == Knots([])
        @test Knots(1:3) == Knots([3, 2, 1])

        k = Knots([4, 5, 6])
        @test ♯(k) == length(k) == 3

        @test Knots([-1, 2, 3]) + 2 * Knots([2, 5]) == Knots([-1, 2, 2, 2, 3, 5, 5])
        Knots([1, 2, 3]) + Knots([2, 4, 5]) == Knots([1, 2, 2, 3, 5])
        2 * Knots([2, 3]) == Knots([2, 2, 3, 3])

        unique(Knots([1, 2, 2, 3])) == Knots([1, 2, 3])

        k = Knots([1, 2, 2, 3])
        @test 𝔫(k, 0.3) == 0
        @test 𝔫(k, 1.0) == 1
        @test 𝔫(k, 2.0) == 2
    end

    @testset "BSplineSpace" begin
        Random.seed!(42)

        P1 = BSplineSpace(2, Knots([1, 3, 5, 6, 8, 9]))
        @test bsplinesupport(2, P1) == 3..8
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
        @test bsplinesupport(i, P) == 12..14

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
        Δ12 = [bsplinebasis(i, P1, t) - sum(A12[i, j] * bsplinebasis(j, P2, t) for j in 1:n2) for i in 1:n1, t in 8 * rand(10)]
        Δ13 = [bsplinebasis(i, P1, t) - sum(A13[i, j] * bsplinebasis(j, P3, t) for j in 1:n3) for i in 1:n1, t in 8 * rand(10)]
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
        Δ45 = [bsplinebasis(i, P4, t) - sum(A45[i, j] * bsplinebasis(j, P5, t) for j in 1:n5) for i in 1:n1, t in 2 * rand(10) .+ 2]
        @test norm(Δ45) < 1e-14
    end

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

    @testset "0th degree basis" begin
        Random.seed!(42)

        p = 0
        k = Knots(rand(10)) + Knots(0, 1)
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in k])

        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in k])
    end

    @testset "1st degree basis" begin
        Random.seed!(42)

        p = 1
        k = Knots(rand(10)) + 2 * Knots(0, 1)
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in k])

        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in k])
    end

    @testset "2nd degree basis" begin
        Random.seed!(42)

        p = 2
        k = Knots(rand(10)) + 3 * Knots(0, 1)
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in k])

        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in k])
    end

    @testset "3rd degree basis" begin
        Random.seed!(42)

        p = 3
        k = Knots(rand(10)) + 4 * Knots(0, 1)
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in k])

        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in k])
    end

    @testset "5th degree basis" begin
        Random.seed!(42)

        p = 5
        k = Knots(rand(10)) + 6 * Knots(0, 1)
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t) for i in 1:n, t in k])

        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) for i in 1:n, t in k])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in rand(5)])
        @test prod([bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t) for i in 1:n, t in k])
    end

    @testset "_bsplinebasis" begin
        Random.seed!(42)

        k = Knots(rand(10).-1) + Knots(rand(10)) + Knots(rand(10).+1)
        ts = rand(10)

        for p in 0:BasicBSpline.MAX_DEGREE
            P = FastBSplineSpace(p,k)
            for t in ts
                j = BasicBSpline._knotindex(P,t)
                _B = BasicBSpline._bsplinebasis(P,t,j)
                B = [bsplinebasis(i,P,t) for i in j-p:j]
                @test norm(collect(_B) - B) < ε
            end
        end
    end

    @testset "BSplineManifold-2dim" begin
        Random.seed!(42)

        P1 = BSplineSpace(1, Knots([0, 0, 1, 1]))
        P2 = BSplineSpace(1, Knots([1, 1, 2, 3, 3]))
        n1 = dim(P1) # 2
        n2 = dim(P2) # 3
        a = [[i, j] for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
        M = BSplineManifold([P1, P2], a)
        @test dim(M) == 2

        P1′ = BSplineSpace(2, Knots([-2, 0, 0, 1, 1, 2]))
        P2′ = BSplineSpace(1, Knots([-3, 1, 2, 1.45, 3, 4]))
        p₊ = [1, 0]
        k₊ = [Knots(), Knots(1.45)]

        @test P1 ⊑ P1′
        @test P2 ⊑ P2′

        M′ = refinement(M, [P1′, P2′])
        M′′ = refinement(M, p₊=p₊, k₊=k₊)
        ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
        @test prod([M(t) ≈ M′(t) for t in ts])
        @test prod([M(t) ≈ M′′(t) for t in ts])
    end

    @testset "FastBSplineManifold-2dim" begin
        Random.seed!(42)

        P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
        P2 = FastBSplineSpace(1, Knots([1, 1, 2, 3, 3]))
        n1 = dim(P1) # 2
        n2 = dim(P2) # 3
        a = [[i, j] for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
        M = FastBSplineManifold([P1, P2], a)
        @test dim(M) == 2

        P1′ = FastBSplineSpace(2, Knots([0, 0, 0, 1, 1, 1]))
        P2′ = FastBSplineSpace(1, Knots([1, 1, 2, 1.45, 3, 3]))
        p₊ = [1, 0]
        k₊ = [Knots(), Knots(1.45)]

        @test P1 ⊆ P1′
        @test P2 ⊆ P2′

        M′ = refinement(M, [P1′, P2′])
        M′′ = refinement(M, p₊=p₊, k₊=k₊)
        ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
        @test prod([M(t) ≈ M′(t) for t in ts])
        @test prod([M(t) ≈ M′′(t) for t in ts])
    end

    @testset "BSplineSurface" begin
        Random.seed!(42)

        P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
        P2 = FastBSplineSpace(1, Knots([1, 1, 2, 3, 3]))
        n1 = dim(P1) # 2
        n2 = dim(P2) # 3
        a = [[i, j] for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
        M = BSplineSurface([P1, P2], a)
        @test dim(M) == 2

        P1′ = FastBSplineSpace(2, Knots([0, 0, 0, 1, 1, 1]))
        P2′ = FastBSplineSpace(1, Knots([1, 1, 2, 1.45, 3, 3]))
        p₊ = [1, 0]
        k₊ = [Knots(), Knots(1.45)]

        @test P1 ⊆ P1′
        @test P2 ⊆ P2′

        M′ = refinement(M, [P1′, P2′])
        M′′ = refinement(M, p₊=p₊, k₊=k₊)
        ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
        @test prod([M(t) ≈ M′(t) for t in ts])
        @test prod([M(t) ≈ M′′(t) for t in ts])
    end

    @testset "Fitting-curve_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + (p1 + 1) * Knots([0, 1])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        a_org = [[i1, rand()] for i1 in 1:n1]
        M = BSplineCurve([P1], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)

        M′ = refinement(M, [P1′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′])
        a_fit = transpose(hcat(a_tmp...))

        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-curve_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + Knots([0, 1]) + p1 * Knots([-rand(), 1 + rand()])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        a_org = [[i1, rand()] for i1 in 1:n1]
        M = BSplineCurve([P1], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1[1+p1:end-p1]) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)

        M′ = refinement(M, [P1′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′])
        a_fit = transpose(hcat(a_tmp...))

        @test norm(a_fit - a_ref) < ε
    end

    @testset "array <-> array of vector" begin
        Random.seed!(42)

        A = rand(3,4)
        A′ = BasicBSpline.array2arrayofvector(A)
        A′′ = BasicBSpline.arrayofvector2array(A′)
        @test A == A′′

        B = rand(3,4,5)
        B′ = BasicBSpline.array2arrayofvector(B)
        B′′ = BasicBSpline.arrayofvector2array(B′)
        @test B == B′′

        C = rand(3,4,5,6)
        C′ = BasicBSpline.array2arrayofvector(C)
        C′′ = BasicBSpline.arrayofvector2array(C′)
        @test C == C′′
    end

    @testset "Fitting-surface_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + (p1 + 1) * Knots([0, 1])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        p2 = 1
        k2 = Knots(rand(4)) + (p2 + 1) * Knots([0, 1])
        P2 = FastBSplineSpace(p2, k2)
        n2 = dim(P2)
        a_org = [[i1, i2, rand()] for i1 in 1:n1, i2 in 1:n2]
        M = BSplineSurface([P1, P2], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)
        p2′ = p2 + 1
        k2′ = k2 + unique(k2) + Knots(rand(2))
        P2′ = FastBSplineSpace(p2′, k2′)

        M′ = refinement(M, [P1′, P2′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′, P2′])
        a_fit = BasicBSpline.arrayofvector2array(a_tmp)

        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-surface_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + Knots([0, 1]) + p1 * Knots([-rand(), 1 + rand()])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        p2 = 1
        k2 = Knots(rand(4)) + Knots([0, 1]) + p2 * Knots([-rand(), 1 + rand()])
        P2 = FastBSplineSpace(p2, k2)
        n2 = dim(P2)
        a_org = [[i1, i2, rand()] for i1 in 1:n1, i2 in 1:n2]
        M = BSplineSurface([P1, P2], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1[1+p1:end-p1]) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)
        p2′ = p2 + 1
        k2′ = k2 + unique(k2[1+p2:end-p2]) + Knots(rand(2))
        P2′ = FastBSplineSpace(p2′, k2′)

        M′ = refinement(M, [P1′, P2′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′, P2′])
        a_fit = BasicBSpline.arrayofvector2array(a_tmp)

        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-solid_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + (p1 + 1) * Knots([0, 1])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        p2 = 1
        k2 = Knots(rand(4)) + (p2 + 1) * Knots([0, 1])
        P2 = FastBSplineSpace(p2, k2)
        n2 = dim(P2)
        p3 = 2
        k3 = Knots(rand(5)) + (p3 + 1) * Knots([0, 1])
        P3 = FastBSplineSpace(p3, k3)
        n3 = dim(P3)
        a_org = [[i1, i2, i3, rand()] for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
        M = BSplineSolid([P1, P2, P3], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)
        p2′ = p2 + 1
        k2′ = k2 + unique(k2) + Knots(rand(2))
        P2′ = FastBSplineSpace(p2′, k2′)
        p3′ = p3 + 1
        k3′ = k3 + unique(k3) + Knots(rand(2))
        P3′ = FastBSplineSpace(p3′, k3′)

        M′ = refinement(M, [P1′, P2′, P3′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′, P2′, P3′])
        a_fit = BasicBSpline.arrayofvector2array(a_tmp)

        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-solid_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = Knots(rand(3)) + Knots([0, 1]) + p1 * Knots([-rand(), 1 + rand()])
        P1 = FastBSplineSpace(p1, k1)
        n1 = dim(P1)
        p2 = 1
        k2 = Knots(rand(4)) + Knots([0, 1]) + p1 * Knots([-rand(), 1 + rand()])
        P2 = FastBSplineSpace(p2, k2)
        n2 = dim(P2)
        p3 = 2
        k3 = Knots(rand(5)) + Knots([0, 1]) + p1 * Knots([-rand(), 1 + rand()])
        P3 = FastBSplineSpace(p3, k3)
        n3 = dim(P3)
        a_org = [[i1, i2, i3, rand()] for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
        M = BSplineSolid([P1, P2, P3], a_org)

        p1′ = p1 + 1
        k1′ = k1 + unique(k1[1+p1:end-p1]) + Knots(rand(2))
        P1′ = FastBSplineSpace(p1′, k1′)
        p2′ = p2 + 1
        k2′ = k2 + unique(k2[1+p2:end-p2]) + Knots(rand(2))
        P2′ = FastBSplineSpace(p2′, k2′)
        p3′ = p3 + 1
        k3′ = k3 + unique(k3[1+p3:end-p3]) + Knots(rand(2))
        P3′ = FastBSplineSpace(p3′, k3′)

        M′ = refinement(M, [P1′, P2′, P3′])
        a_ref = M′.controlpoints

        a_tmp = fittingcontrolpoints(M, [P1′, P2′, P3′])
        a_fit = BasicBSpline.arrayofvector2array(a_tmp)

        @test norm(a_fit - a_ref) < 1.0e-6
    end

end
