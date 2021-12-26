@testset "Derivative" begin
    p_max = 5
    ts = rand(10)
    dt = 1e-7

    for p in 1:p_max
        k = Knots(rand(20)) + (p+1)*Knots(0,1)
        P = BSplineSpace{p}(k)
        P0 = BSplineDerivativeSpace{0}(P)
        for t in ts, i in 1:dim(P)
            @test bsplinebasis(P,i,t) == bsplinebasis(P0,i,t)
        end
        
        for r in 1:p_max
            P1 = BSplineDerivativeSpace{r-1}(P)
            P2 = BSplineDerivativeSpace{r}(P)        
            for t in ts, i in 1:dim(P)
                d1 = (bsplinebasis(P1,6,t+dt) - bsplinebasis(P1,6,t-dt))/2dt
                d2 = bsplinebasis(P2,6,t)
                @test d1 ≈ d2 atol=1e-4
            end
        end
    end
end
