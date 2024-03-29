@testset "Plots" begin
    dir_out = joinpath(@__DIR__, "out_plot")
    rm(dir_out; force=true, recursive=true)
    mkpath(dir_out)

    @testset "B-spline space" begin
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        pl = plot(P)
        path_img = joinpath(dir_out, "bspline_space.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline space with knotvector" begin
        k = knotvector"112111113111111121114"
        P = BSplineSpace{3}(k)
        pl = plot(P; label="B-spline basis"); plot!(k; label="knot vector")
        path_img = joinpath(dir_out, "bspline_space_with_knotvector.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "Uniform B-spline space" begin
        p = 3
        k = UniformKnotVector(0:8)
        P = BSplineSpace{p}(k)
        pl = plot(P)
        path_img = joinpath(dir_out, "bspline_space_uniform.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "Derivative B-spline space" begin
        p = 3
        k = UniformKnotVector(0:8)
        dP = BSplineDerivativeSpace{1}(BSplineSpace{p}(k))
        pl = plot(dP)
        path_img = joinpath(dir_out, "bspline_space_derivative.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline spaces with knotvectors" begin
        k1 = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
        P1 = BSplineSpace{3}(k1)
        P2 = expandspace_R(P1, KnotVector([2.0, 3.0, 5.5]))

        pl = plot(P1, label="P1", color=:red)
        plot!(P2, label="P2", color=:green)
        plot!(knotvector(P1), label="k1", gap=0.01, color=:red)
        plot!(knotvector(P2); label="k2", gap=0.01, offset=-0.03, color=:green)

        path_img = joinpath(dir_out, "bspline_spaces_and_their_knotvectors.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "Tensor product of B-spline spaces" begin
        # Define piecewise polynomial spaces
        k1 = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
        k2 = knotvector"31 2  121"
        P1 = BSplineSpace{3}(k1)
        P2 = BSplineSpace{2}(k2)

        # Choose indices
        i1 = 3
        i2 = 4

        # Visualize basis functions
        xs = range(0,10,length=100)
        ys = range(1,9,length=100)
        pl = surface(xs, ys, bsplinebasis.(P1,i1,xs') .* bsplinebasis.(P2,i2,ys))
        plot!(P1; plane=:xz, label="P1", color=:red)
        plot!(P2; plane=:yz, label="P2", color=:green)
        plot!(k1; plane=:xz, label="k1", color=:red, markersize=2)
        plot!(k2; plane=:yz, label="k2", color=:green, markersize=2)

        path_img = joinpath(dir_out, "bspline_space_tensor_product.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, P)
        pl = plot(M)
        path_img = joinpath(dir_out, "bspline_curve_2d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = BSplineManifold(a, P)
        pl = plot(M)
        path_img = joinpath(dir_out, "bspline_curve_3d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline surface in 2d" begin
        k1 = KnotVector(1:8)
        k2 = KnotVector(1:10)
        P1 = BSplineSpace{2}(k1)
        P2 = BSplineSpace{2}(k2)
        n1 = dim(P1)
        n2 = dim(P2)
        a = [SVector(i+j/2-cos(i-j)/5, j-i/2+sin(i+j)/5) for i in 1:n1, j in 1:n2]
        M = BSplineManifold(a,(P1,P2))
        pl = plot(M)
        path_img = joinpath(dir_out, "bspline_surface_2d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
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
        path_img = joinpath(dir_out, "bspline_surface_3d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "B-spline solid in 3d" begin
        P1 = P2 = BSplineSpace{3}(KnotVector([0,0,0,0,1,1,1,1]))
        P3 = BSplineSpace{3}(KnotVector([0,0,0,0,1,2,3,4,5,6,6,6,6]))
        n1 = dim(P1)
        n2 = dim(P2)
        n3 = dim(P3)
        a = [SVector(i+j/2-cos(i-j)/5+sin(k), j-i/2+sin(i+j)/5+cos(k), k) for i in 1:n1, j in 1:n2, k in 1:n3]
        M = BSplineManifold(a, P1, P2, P3);
        pl = plot(M; controlpoints=(;markersize=1))
        path_img = joinpath(dir_out, "bspline_solid_3d.png")
        @test !isfile(path_img)
        # https://github.com/hyrodium/BasicBSpline.jl/pull/374#issuecomment-1927050884
        # savefig(pl, path_img)
        # @test isfile(path_img)
    end

    @testset "Rational B-spline curve in 2d" begin
        a = [SVector(0, 0), SVector(1, 1), SVector(2, -1), SVector(3, 0), SVector(4, -2), SVector(5, 1)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, P)
        pl = plot(M)
        path_img = joinpath(dir_out, "rational_bspline_curve_2d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "Rational B-spline curve in 3d" begin
        a = [SVector(0, 0, 2), SVector(1, 1, 1), SVector(2, -1, 4), SVector(3, 0, 0), SVector(4, -2, 0), SVector(5, 1, 2)]
        w = rand(6)
        p = 3
        k = KnotVector(0:3)+p*KnotVector([0,3])
        P = BSplineSpace{p}(k)
        M = RationalBSplineManifold(a, w, P)
        pl = plot(M)
        path_img = joinpath(dir_out, "rational_bspline_curve_3d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
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
        path_img = joinpath(dir_out, "rational_bspline_surface_3d.png")
        @test !isfile(path_img)
        # https://github.com/hyrodium/BasicBSpline.jl/pull/374#issuecomment-1927050884
        # savefig(pl, path_img)
        # @test isfile(path_img)
    end
end
