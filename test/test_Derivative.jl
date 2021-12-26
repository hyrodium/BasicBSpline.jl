@testset "Derivative" begin
    Random.seed!(42)
    p_max = 4
    ts = rand(10)
    dt = 1e-7

    # Not sure why this @testset doesn't work fine.
    # @testset "$(p)-th degree basis" for p in 0:p_max
    for p in 0:p_max
        k = Knots(rand(20)) + (p+1)*Knots(0,1)
        P = BSplineSpace{p}(k)
        P0 = BSplineDerivativeSpace{0}(P)
        for t in ts, i in 1:dim(P)
            @test bsplinebasis(P,i,t) == bsplinebasis(P0,i,t)
            @test bsplinebasis₋₀(P,i,t) == bsplinebasis₋₀(P0,i,t)
            @test bsplinebasis₊₀(P,i,t) == bsplinebasis₊₀(P0,i,t)
        end
        @test dim(P) == dim(P0)

        # Check the values of the derivative of B-spline basis funcitons
        n = dim(P0)
        for t in k
            s = sum([bsplinebasis(P0, i, t) for i in 1:n])
            s₊₀ = sum([bsplinebasis₊₀(P0, i, t) for i in 1:n])
            s₋₀ = sum([bsplinebasis₋₀(P0, i, t) for i in 1:n])
            if t == k[1]
                @test s ≈ 1
                @test s₊₀ ≈ 1
                @test s₋₀ == 0
            elseif t == k[end]
                @test s ≈ 1
                @test s₊₀ == 0
                @test s₋₀ ≈ 1
            else
                @test s ≈ 1
                @test s₊₀ ≈ 1
                @test s₋₀ ≈ 1
            end
        end

        for r in 1:p_max
            P1 = BSplineDerivativeSpace{r-1}(P)
            P2 = BSplineDerivativeSpace{r}(P)
            @test degree(P2) == p-r
            @test dim(P2) == dim(P)-r
            for t in ts, i in 1:dim(P)
                d1 = (bsplinebasis(P1,i,t+dt) - bsplinebasis(P1,i,t-dt))/2dt
                d2 = bsplinebasis(P2,i,t)
                @test d1 ≈ d2
            end
        end
    end
end
