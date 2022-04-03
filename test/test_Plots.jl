@testset "Plots" begin
    @testset "B-spline space" begin
        p = 3
        k = KnotVector(0:3)+p*KnotVector(0,3)
        P = BSplineSpace{p}(k)
        plot(P)
    end

    @testset "uniform B-spline space" begin
        p = 3
        k = UniformKnotVector(0:8)
        P = UniformBSplineSpace{p}(k)
        plot(P)
    end

    @testset "B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector(0,3)
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, (P,))
        plot(M)
    end

    @testset "B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector(0,3)
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, (P,))
        plot(M)
    end

    @testset "Rational B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector(0,3)
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, (P,))
        plot(M)
    end

    @testset "Rational B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector(0,3)
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, (P,))
        plot(M)
    end
end
