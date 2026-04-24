@testset "PeriodicBSplineBasis" begin
    integer_knotsets = (
        KnotVector([0, 1, 2, 3, 4]),
        UniformKnotVector(0:6),
    )
    float_knotsets = (
        KnotVector([0.0, 1.0, 2.0, 3.0, 4.0]),
        KnotVector([0.0, 0.3, 0.6, 1.0]),
        KnotVector([-1.0, 0.25, 0.75, 1.5, 2.0, 3.1]),
    )

    # Shared checks across a range of sample points. Generic over knotset element
    # type; the caller decides which degrees are appropriate.
    function _periodic_basis_checks(P, knot_floats)
        n = dim(P)
        L = BasicBSpline.period(P)
        kfirst, klast = knot_floats[1], knot_floats[end]

        # Periodicity: B_i(t) ≈ B_i(t ± L)
        for t in range(kfirst - 0.3, klast + 0.3, length=13)
            for i in 1:n
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(P, i, t + L) atol = 1.0e-15
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(P, i, t - L) atol = 1.0e-15
            end
        end

        # Partition of unity across one period
        for t in range(kfirst, klast; length=17)[1:end-1]
            @test sum(bsplinebasis(P, i, t) for i in 1:n) ≈ 1
        end

        # Non-negativity across an extended range
        for t in range(kfirst - 1, klast + 1, length=21)
            for i in 1:n
                @test bsplinebasis(P, i, t) ≥ -1e-12
            end
        end

        # intervalindex always lands in 1:n
        for t in range(kfirst - 2, klast + 2, length=9)
            j = intervalindex(P, t)
            @test 1 ≤ j ≤ n
        end

        # bsplinebasis == bsplinebasis₊₀ on the circle (no boundary to worry about)
        for t in range(kfirst + 0.01, klast - 0.01, length=7)
            for i in 1:n
                @test bsplinebasis(P, i, t) == bsplinebasis₊₀(P, i, t)
            end
        end
    end

    @testset "integer knots (p = $p)" for p in 0:3
        for k in integer_knotsets
            length(k) < p + 2 && continue
            P = PeriodicBSplineSpace{p}(k)
            _periodic_basis_checks(P, [float(k[j]) for j in 1:length(k)])
        end
    end

    @testset "float knots (p = $p)" for p in 1:3
        for k in float_knotsets
            length(k) < p + 2 && continue
            P = PeriodicBSplineSpace{p}(k)
            _periodic_basis_checks(P, [float(k[j]) for j in 1:length(k)])
        end
    end

    @testset "wrap-around basis values" begin
        # Uniform periodic quadratic: basis B_4 must equal a shifted copy of B_1.
        k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
        P = PeriodicBSplineSpace{2}(k)
        for t in range(-2.0, 6.0, length=17)
            @test bsplinebasis(P, 4, t) ≈ bsplinebasis(P, 1, t - 3)
        end
    end

    @testset "support" begin
        k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
        P = PeriodicBSplineSpace{2}(k)
        @test bsplinesupport_R(P, 1) == 0.0 .. 3.0
        @test bsplinesupport_R(P, 2) == 1.0 .. 4.0
        @test bsplinesupport_R(P, 4) == 3.0 .. 6.0   # wraps past the period
        @test bsplinesupport_I(P, 4) == bsplinesupport_R(P, 4)
    end
end
