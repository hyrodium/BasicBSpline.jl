@testset "ChainRules" begin
    k = KnotVector(rand(20))
    p = 3
    P = BSplineSpace{p}(k)
    dP0 = BSplineDerivativeSpace{0}(P)
    dP1 = BSplineDerivativeSpace{1}(P)
    dP2 = BSplineDerivativeSpace{2}(P)
    @testset "bsplinebasis" begin
        for _ in 1:10
            t = rand(domain(P))
            for _P in (P, dP0, dP1, dP2), i in 1:dim(_P)
                test_frule(bsplinebasis, _P, i, t)
                test_rrule(bsplinebasis, _P, i, t)
                test_frule(bsplinebasis₊₀, _P, i, t)
                test_rrule(bsplinebasis₊₀, _P, i, t)
                test_frule(bsplinebasis₋₀, _P, i, t)
                test_rrule(bsplinebasis₋₀, _P, i, t)
            end
        end
    end
    @testset "bsplinebasisall" begin
        for _ in 1:10
            t = rand(domain(_P))
            for _P in (P, dP0, dP1, dP2), i in 1:length(k)-2p-1
                test_frule(bsplinebasisall, _P, i, t)
                test_rrule(bsplinebasisall, _P, i, t)
            end
        end
    end
end
