@testset "Fitting" begin
    ε = 1.0e-8

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

        a_tmp = fittingcontrolpoints(M, [P1′], domain=:R)
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

        a_tmp = fittingcontrolpoints(M, [P1′, P2′], domain=:R)
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

        a_tmp = fittingcontrolpoints(M, [P1′, P2′, P3′], domain=:R)
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
