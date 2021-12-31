@testset "Derivative" begin
    Random.seed!(42)
    p_max = 4
    ts = rand(10)
    dt = 1e-8

    # Not sure why this @testset doesn't work fine.
    # @testset "$(p)-th degree basis" for p in 0:p_max
    for p in 0:p_max
        k = KnotVector(rand(20)) + (p+1)*KnotVector(0,1)
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

        P1 = BSplineDerivativeSpace{1}(P)
        for t in ts, i in 1:dim(P)
            @test bsplinebasis′(P,i,t) == bsplinebasis(P1,i,t)
        end

        for r in 1:p_max
            Pa = BSplineDerivativeSpace{r}(P)
            Pb = BSplineDerivativeSpace{r-1}(P)
            @test degree(Pa) == p-r
            @test dim(Pa) == dim(P)
            @test properdim(Pa) == properdim(P)-r
            @test domain(P) == domain(Pa) == domain(Pb)
            for t in ts, i in 1:dim(P)
                d1 = bsplinebasis(Pa,i,t)
                d2 = (bsplinebasis(Pb,i,t+dt) - bsplinebasis(Pb,i,t-dt))/2dt
                d3 = bsplinebasis′(Pb,i,t)
                @test d1 ≈ d2 rtol=1e-7
                @test d1 == d3
            end
        end
    end

    @testset "bsplinebasisall" begin
        Random.seed!(42)

        k = KnotVector(rand(10).-1) + KnotVector(rand(10)) + KnotVector(rand(10).+1)
        ts = rand(10)

        for p in 0:5
            P = BSplineSpace{p}(k)

            for r in 0:5
                dP = BSplineDerivativeSpace{r}(P)

                for t in ts
                    j = intervalindex(dP,t)
                    B = collect(bsplinebasisall(dP,j,t))

                    _B = bsplinebasis.(dP,j:j+p,t)
                    @test _B ≈ B

                    _B = bsplinebasis₊₀.(dP,j:j+p,t)
                    @test _B ≈ B

                    _B = bsplinebasis₋₀.(dP,j:j+p,t)
                    @test _B ≈ B
                end
            end
        end
    end

    @testset "Rational" begin
        k = KnotVector{Int}(1:12)
        P = BSplineSpace{3}(k)
        dP = BSplineDerivativeSpace{1}(P)
        k isa KnotVector{Int}
        dP isa BSplineDerivativeSpace{1,BSplineSpace{3,Int}}
        bsplinebasis(dP,1,11//5) isa Rational{Int}
        bsplinebasis₊₀(dP,1,11//5) isa Rational{Int}
        bsplinebasis₋₀(dP,1,11//5) isa Rational{Int}

        bsplinebasis(dP,1,11//5)   ===
        bsplinebasis₊₀(dP,1,11//5) ===
        bsplinebasis₋₀(dP,1,11//5) ===
        bsplinebasis′(P,1,11//5)   ===
        bsplinebasis′₊₀(P,1,11//5) ===
        bsplinebasis′₋₀(P,1,11//5) === 16//25
    end

    @testset "Endpoints" begin
        p = 2

        k = KnotVector(1:2) + (p+1)*KnotVector(0,3)
        P0 = BSplineSpace{0}(k)
        P1 = BSplineSpace{1}(k)
        P2 = BSplineSpace{2}(k)

        @test !isproper(P0)
        @test !isproper(P1)
        @test isproper(P2)

        n0 = dim(P0)
        n1 = dim(P1)
        n2 = dim(P2)

        @test bsplinebasis′.(P0,1:n0,0) == bsplinebasis′₊₀.(P0,1:n0,0) == [0,0,0,0,0,0,0]
        @test bsplinebasis′.(P1,1:n1,0) == bsplinebasis′₊₀.(P1,1:n1,0) == [0,-1,1,0,0,0]
        @test bsplinebasis′.(P2,1:n2,0) == bsplinebasis′₊₀.(P2,1:n2,0) == [-2,2,0,0,0]
        @test bsplinebasis′.(P0,1:n0,3) == bsplinebasis′₋₀.(P0,1:n0,3) == [0,0,0,0,0,0,0]
        @test bsplinebasis′.(P1,1:n1,3) == bsplinebasis′₋₀.(P1,1:n1,3) == [0,0,0,-1,1,0]
        @test bsplinebasis′.(P2,1:n2,3) == bsplinebasis′₋₀.(P2,1:n2,3) == [0,0,0,-2,2]

        @test bsplinebasis′₋₀.(P0,1:n0,0) == [0,0,0,0,0,0,0]
        @test bsplinebasis′₋₀.(P1,1:n1,0) == [0,0,0,0,0,0]
        @test bsplinebasis′₋₀.(P2,1:n2,0) == [0,0,0,0,0]
        @test bsplinebasis′₊₀.(P0,1:n0,3) == [0,0,0,0,0,0,0]
        @test bsplinebasis′₊₀.(P1,1:n1,3) == [0,0,0,0,0,0]
        @test bsplinebasis′₊₀.(P2,1:n2,3) == [0,0,0,0,0]
    end
end
