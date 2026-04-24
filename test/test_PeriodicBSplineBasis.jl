@testset "PeriodicBSplineBasis" begin
    knotsets = (
        KnotVector([0.0, 1.0, 2.0, 3.0, 4.0]),
        KnotVector([0.0, 0.3, 0.6, 1.0]),
        KnotVector([-1.0, 0.25, 0.75, 1.5, 2.0, 3.1]),
        UniformKnotVector(0:6),
    )
    degrees = 0:3

    @testset "degree p = $p" for p in degrees
        for k in knotsets
            length(k) < p + 2 && continue
            P = PeriodicBSplineSpace{p}(k)
            n = dim(P)
            L = BasicBSpline.period(P)

            # Periodicity: B_i(t) == B_i(t + L) for any t and i
            for t in range(float(k[1]) - 0.3, float(k[end]) + 0.3, length=13)
                for i in 1:n
                    @test bsplinebasis(P, i, t) ≈ bsplinebasis(P, i, t + L) atol = 1.0e-15
                    @test bsplinebasis(P, i, t) ≈ bsplinebasis(P, i, t - L) atol = 1.0e-15
                end
            end

            # Partition of unity across one period
            for t in range(float(k[1]), float(k[end]); length=17)[1:end-1]
                s = sum(bsplinebasis(P, i, t) for i in 1:n)
                @test s ≈ 1
            end

            # Non-negativity
            for t in range(float(k[1]) - 1, float(k[end]) + 1, length=21)
                for i in 1:n
                    @test bsplinebasis(P, i, t) ≥ -1e-12
                end
            end

            # intervalindex is in 1:n
            for t in range(float(k[1]) - 2, float(k[end]) + 2, length=9)
                j = intervalindex(P, t)
                @test 1 ≤ j ≤ n
            end

            # bsplinebasis == bsplinebasis₊₀ on the circle
            for t in range(float(k[1]) + 0.01, float(k[end]) - 0.01, length=7)
                for i in 1:n
                    @test bsplinebasis(P, i, t) == bsplinebasis₊₀(P, i, t)
                end
            end
        end
    end

    @testset "wrap-around basis values" begin
        # Uniform periodic quadratic with knots [0,1,2,3,4], period 4, n=4.
        # Basis B_4 has support cyclically [3,6] = [3,4] ∪ [0,2].
        # By uniform-periodic symmetry, B_4(t) on [0,2] ∪ [3,4] mirrors B_1 on a shift.
        k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
        P = PeriodicBSplineSpace{2}(k)
        for t in range(-2.0, 6.0, length=17)
            # B_4(t) must equal B_1(t - 3) by a period-preserving shift argument.
            @test bsplinebasis(P, 4, t) ≈ bsplinebasis(P, 1, t - 3)
        end
    end

    @testset "support" begin
        k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
        P = PeriodicBSplineSpace{2}(k)
        @test bsplinesupport_R(P, 1) == 0.0 .. 3.0
        @test bsplinesupport_R(P, 2) == 1.0 .. 4.0
        @test bsplinesupport_R(P, 4) == 3.0 .. 6.0   # wraps
        @test bsplinesupport_I(P, 4) == bsplinesupport_R(P, 4)
    end
end
