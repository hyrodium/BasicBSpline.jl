@testset "ChainRules" begin
    k = KnotVector(rand(12))
    p = 2
    P = BSplineSpace{p}(k)
    @testset "bsplinebasis" begin
        for _ in 1:10
            t = rand(domain(P))
            for i in 1:dim(P)
                test_frule(bsplinebasis, P, i, t)
                test_rrule(bsplinebasis, P, i, t)
                test_frule(bsplinebasis₊₀, P, i, t)
                test_rrule(bsplinebasis₊₀, P, i, t)
                test_frule(bsplinebasis₋₀, P, i, t)
                test_rrule(bsplinebasis₋₀, P, i, t)
            end
            for i in 1:length(k)-2p-1
                test_frule(bsplinebasisall, P, i, t)
                test_rrule(bsplinebasisall, P, i, t)
            end
        end
    end
end
