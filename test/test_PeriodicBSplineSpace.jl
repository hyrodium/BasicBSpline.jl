@testset "PeriodicBSplineSpace" begin
    @testset "constructor" begin
        P1 = PeriodicBSplineSpace{2}(KnotVector(1:5))
        P2 = PeriodicBSplineSpace{2}(P1)
        P3 = PeriodicBSplineSpace{2,Int}(P1)
        P4 = PeriodicBSplineSpace{2,Real}(P1)
        P5 = PeriodicBSplineSpace{2,Float64}(P1)
        P6 = PeriodicBSplineSpace{2,Float64}(KnotVector(1:5))
        P7 = PeriodicBSplineSpace{2,Float64,KnotVector{Float64}}(KnotVector(1:5))
        @test P1 == P2 == P3 == P4 == P5 == bsplinespace(P1) == bsplinespace(P2)
        @test P5 == P6 == P7
        @test typeof(P5) === typeof(P6) === typeof(P7)
        @test P1 isa PeriodicBSplineSpace{2,Int}
        @test P5 isa PeriodicBSplineSpace{2,Float64}
        @test P1 === PeriodicBSplineSpace(P1)
        @test_throws DomainError PeriodicBSplineSpace{-1}(KnotVector(1:5))
        # Too few knots for the degree (require length(k) ≥ p+2)
        @test_throws DomainError PeriodicBSplineSpace{3}(KnotVector(1:3))
        # AbstractBSplineSpace hierarchy
        @test P1 isa BasicBSpline.AbstractBSplineSpace{2,Int}
    end

    @testset "dim, domain, period" begin
        P = PeriodicBSplineSpace{2}(KnotVector([0.0, 1.0, 2.0, 3.0, 4.0]))
        @test dim(P) == 4                       # length(k) - 1, independent of p
        @test domain(P) == 0.0 .. 4.0
        @test BasicBSpline.period(P) ≈ 4.0
        @test degree(P) == 2

        # dim is independent of p
        for p in 0:4
            Q = PeriodicBSplineSpace{p}(KnotVector(0:6))
            @test dim(Q) == 6
            @test BasicBSpline.period(Q) == 6
        end
    end

    @testset "equality, copy, hash" begin
        k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
        P1 = PeriodicBSplineSpace{2}(k)
        P2 = PeriodicBSplineSpace{2}(copy(k))
        P3 = PeriodicBSplineSpace{3}(k)
        @test P1 == P2
        @test hash(P1) == hash(P2)
        @test P1 != P3
        Pc = copy(P1)
        @test Pc == P1
        @test Pc !== P1
    end

    @testset "UniformKnotVector" begin
        kU = UniformKnotVector(0:4)
        P = PeriodicBSplineSpace{2}(kU)
        @test dim(P) == 4
        @test BasicBSpline.period(P) == 4
        @test P == PeriodicBSplineSpace{2}(KnotVector(0:4))
    end

    @testset "_periodic_knot / _periodic_reduce" begin
        k = KnotVector([0.0, 1.0, 2.5, 3.0, 5.0])
        P = PeriodicBSplineSpace{2}(k)
        L = BasicBSpline.period(P)
        # Identity on 1:n
        for j in 1:4
            @test BasicBSpline._periodic_knot(P, j) == k[j]
        end
        # Cyclic extension
        @test BasicBSpline._periodic_knot(P, 5) == k[1] + L
        @test BasicBSpline._periodic_knot(P, 6) == k[2] + L
        @test BasicBSpline._periodic_knot(P, 0) == k[4] - L
        @test BasicBSpline._periodic_knot(P, -1) == k[3] - L
        # _periodic_reduce
        @test BasicBSpline._periodic_reduce(P, 0.5) == 0.5
        @test BasicBSpline._periodic_reduce(P, 0.5 + L) ≈ 0.5
        @test BasicBSpline._periodic_reduce(P, 0.5 - L) ≈ 0.5
    end
end
