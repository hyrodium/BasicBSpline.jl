@testset "BSplineBasis" begin
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
end
