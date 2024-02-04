@testset "Plots" begin
    dir_out = joinpath(@__DIR__, "out_plot")
    rm(dir_out; force=true, recursive=true)
    mkpath(dir_out)

    @testset "B-spline space" begin
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        pl = plot(P)
        savefig(pl, joinpath(dir_out, "bspline_space.png"))
    end

    @testset "Uniform B-spline space" begin
        p = 3
        k = UniformKnotVector(0:8)
        P = BSplineSpace{p}(k)
        pl = plot(P)
        savefig(pl, joinpath(dir_out, "bspline_space_uniform.png"))
    end

    @testset "Derivative B-spline space" begin
        p = 3
        k = UniformKnotVector(0:8)
        dP = BSplineDerivativeSpace{1}(BSplineSpace{p}(k))
        pl = plot(dP)
        savefig(pl, joinpath(dir_out, "bspline_space_derivative.png"))
    end

    @testset "B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, P)
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "bspline_curve_2d.png"))
    end

    @testset "B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, P)
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "bspline_curve_3d.png"))
    end

    @testset "B-spline surface in 3d" begin
        k1 = KnotVector(1:8)
        k2 = KnotVector(1:10)
        P1 = BSplineSpace{2}(k1)
        P2 = BSplineSpace{2}(k2)
        n1 = dim(P1)
        n2 = dim(P2)
        a = [SVector(i,j,sin(i+j)) for i in 1:n1, j in 1:n2]
        M = BSplineManifold(a,(P1,P2))
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "bspline_surface_3d.png"))
    end

    @testset "Rational B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, P)
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "rational_bspline_curve_2d.png"))
    end

    @testset "Rational B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, P)
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "rational_bspline_curve_3d.png"))
    end

    @testset "Rational B-spline surface in 3d" begin
        k1 = KnotVector(1:8)
        k2 = KnotVector(1:10)
        P1 = BSplineSpace{2}(k1)
        P2 = BSplineSpace{2}(k2)
        n1 = dim(P1)
        n2 = dim(P2)
        a = [SVector(i,j,sin(i+j)) for i in 1:n1, j in 1:n2]
        w = rand(n1,n2)
        M = RationalBSplineManifold(a,w,(P1,P2))
        pl = plot(M)
        savefig(pl, joinpath(dir_out, "rational_bspline_surface_3d.png"))
    end
end
