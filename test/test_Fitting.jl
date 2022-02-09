@testset "Fitting" begin
    ε = 1.0e-8

    @testset "Fitting_tinycase" begin
        Random.seed!(42)

        p = 2
        k = KnotVector(1:2p+2)
        P = BSplineSpace{p}(k)
        n = dim(P)
        @test n == p+1
        @test IntervalSets.width(domain(P)) == 1

        a_org = [rand() for i in 1:n]
        f(t) = sum(bsplinebasis(P,i,t)*a_org[i] for i in 1:n)

        a_fit = fittingcontrolpoints(f, (P,), domain=:I)
        @test norm(a_fit - a_org) < ε

        a_fit = fittingcontrolpoints(f, (P,), domain=:R)
        @test norm(a_fit - a_org) < ε
    end

    @testset "Fitting-curve_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + (p1 + 1) * KnotVector([0, 1])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        a_org = [Point(i1, rand()) for i1 in 1:n1]
        M = BSplineManifold(a_org, (P1,))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′,))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′,))
        @test norm(a_fit - a_ref) < ε

        a_fit = fittingcontrolpoints(M, (P1′,), domain=:R)
        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-curve_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + KnotVector([0, 1]) + p1 * KnotVector([-rand(), 1 + rand()])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        a_org = [Point(i1, rand()) for i1 in 1:n1]
        M = BSplineManifold(a_org, (P1,))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′,))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′,))

        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-surface_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + (p1 + 1) * KnotVector([0, 1])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        p2 = 1
        k2 = KnotVector(rand(4)) + (p2 + 1) * KnotVector([0, 1])
        P2 = BSplineSpace{p2}(k2)
        n2 = dim(P2)
        a_org = [Point(i1, i2, rand()) for i1 in 1:n1, i2 in 1:n2]
        M = BSplineManifold(a_org, (P1, P2))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))
        P2′ = expandspace(P2, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′, P2′))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′, P2′))
        @test norm(a_fit - a_ref) < ε

        a_fit = fittingcontrolpoints(M, (P1′, P2′), domain=:R)
        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-surface_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + KnotVector([0, 1]) + p1 * KnotVector([-rand(), 1 + rand()])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        p2 = 1
        k2 = KnotVector(rand(4)) + KnotVector([0, 1]) + p2 * KnotVector([-rand(), 1 + rand()])
        P2 = BSplineSpace{p2}(k2)
        n2 = dim(P2)
        a_org = [Point(i1, i2, rand()) for i1 in 1:n1, i2 in 1:n2]
        M = BSplineManifold(a_org, (P1, P2))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))
        P2′ = expandspace(P2, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′, P2′))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′, P2′))
        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-solid_R" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + (p1 + 1) * KnotVector([0, 1])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        p2 = 1
        k2 = KnotVector(rand(4)) + (p2 + 1) * KnotVector([0, 1])
        P2 = BSplineSpace{p2}(k2)
        n2 = dim(P2)
        p3 = 2
        k3 = KnotVector(rand(5)) + (p3 + 1) * KnotVector([0, 1])
        P3 = BSplineSpace{p3}(k3)
        n3 = dim(P3)
        a_org = [Point(i1, i2, i3, rand()) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
        M = BSplineManifold(a_org, (P1, P2, P3))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))
        P2′ = expandspace(P2, p₊=1, k₊=KnotVector(rand(2)))
        P3′ = expandspace(P3, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′, P2′, P3′))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′, P2′, P3′))
        @test norm(a_fit - a_ref) < ε

        a_fit = fittingcontrolpoints(M, (P1′, P2′, P3′), domain=:R)
        @test norm(a_fit - a_ref) < ε
    end

    @testset "Fitting-solid_I" begin
        Random.seed!(42)

        p1 = 2
        k1 = KnotVector(rand(3)) + KnotVector([0, 1]) + p1 * KnotVector([-rand(), 1 + rand()])
        P1 = BSplineSpace{p1}(k1)
        n1 = dim(P1)
        p2 = 1
        k2 = KnotVector(rand(4)) + KnotVector([0, 1]) + p1 * KnotVector([-rand(), 1 + rand()])
        P2 = BSplineSpace{p2}(k2)
        n2 = dim(P2)
        p3 = 2
        k3 = KnotVector(rand(5)) + KnotVector([0, 1]) + p1 * KnotVector([-rand(), 1 + rand()])
        P3 = BSplineSpace{p3}(k3)
        n3 = dim(P3)
        a_org = [Point(i1, i2, i3, rand()) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
        M = BSplineManifold(a_org, (P1, P2, P3))

        P1′ = expandspace(P1, p₊=1, k₊=KnotVector(rand(2)))
        P2′ = expandspace(P2, p₊=1, k₊=KnotVector(rand(2)))
        P3′ = expandspace(P3, p₊=1, k₊=KnotVector(rand(2)))

        M′ = refinement(M, (P1′, P2′, P3′))
        a_ref = controlpoints(M′)

        a_fit = fittingcontrolpoints(M, (P1′, P2′, P3′))
        @test norm(a_fit - a_ref) < 1.0e-6
    end

    @testset "uniform" begin
        k = UniformKnotVector(1:42)
        for p in 0:3
            P = UniformBSplineSpace{p}(k)
            @test BasicBSpline.innerproduct_R(P) ≈ BasicBSpline.innerproduct_R(BSplineSpace(P))
            @test BasicBSpline.innerproduct_R(P) isa Symmetric{Float64,Matrix{Float64}}
            @test BasicBSpline.innerproduct_R(BSplineSpace(P)) isa Symmetric{Float64,Matrix{Float64}}
        end

        k = UniformKnotVector(1:2:42)
        for p in 0:3
            P = UniformBSplineSpace{p}(k)
            @test BasicBSpline.innerproduct_R(P) ≈ BasicBSpline.innerproduct_R(BSplineSpace(P))
            @test BasicBSpline.innerproduct_R(P) isa Symmetric{Float64,Matrix{Float64}}
            @test BasicBSpline.innerproduct_R(BSplineSpace(P)) isa Symmetric{Float64,Matrix{Float64}}
        end

        k = UniformKnotVector(1:2:42//1)
        for p in 0:3
            P = UniformBSplineSpace{p}(k)
            @test BasicBSpline.innerproduct_R(P) ≈ BasicBSpline.innerproduct_R(BSplineSpace(P))
            @test BasicBSpline.innerproduct_R(P) isa Symmetric{Rational{Int},Matrix{Rational{Int}}}
            @test BasicBSpline.innerproduct_I(P) isa Symmetric{Float64,Matrix{Float64}}  # TODO: This should be Rational
            @test BasicBSpline.innerproduct_R(BSplineSpace(P)) isa Symmetric{Float64,Matrix{Float64}}
        end
    end
end
