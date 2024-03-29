@testset "Derivative" begin
    p_max = 4
    ts = rand(10)
    dt = 1e-8

    @testset "constructor" begin
        P1 = BSplineSpace{2}(KnotVector(1:8))
        P2 = BSplineSpace{2}(UniformKnotVector(1:8))
        P3 = BSplineSpace{2}(UniformKnotVector(1:1:8))
        dP1 = BSplineDerivativeSpace{1}(P1)
        dP2 = BSplineDerivativeSpace{1}(P2)
        dP3 = BSplineDerivativeSpace{1}(P3)
        dP4 = BSplineDerivativeSpace{2}(P3)
        dP5 = BSplineDerivativeSpace{2,typeof(P2)}(P2)
        dP6 = BSplineDerivativeSpace{2,typeof(P3)}(P2)
        dP7 = BSplineDerivativeSpace{2,typeof(P3)}(dP6)
        @test dP1 == dP2 == dP3 != dP4
        @test dP1 !== dP2
        @test hash(dP1) == hash(dP2) == hash(dP3) != hash(dP4)
        @test_throws MethodError BSplineDerivativeSpace{1,typeof(P2)}(dP1)
        @test_throws MethodError BSplineDerivativeSpace(dP1)
        @test_throws MethodError BSplineDerivativeSpace{}(dP1)
        @test dP3 === BSplineDerivativeSpace{1,typeof(P3)}(dP2)
        @test dP2 !== BSplineDerivativeSpace{1,typeof(P3)}(dP2)
        @test dP2 === BSplineDerivativeSpace{1,typeof(P2)}(dP2)
        @test dP1 isa BSplineDerivativeSpace{1,<:BSplineSpace{p,T,KnotVector{T}} where {p,T}}
        @test dP2 isa BSplineDerivativeSpace{1,<:UniformBSplineSpace{p,T} where {p,T}}
        @test dP4 == dP5 == dP6 == dP7
        @test hash(dP4) == hash(dP5) == hash(dP6) == hash(dP7)

        @test dP2 == convert(BSplineDerivativeSpace,dP2)
        @test dP2 == convert(BSplineDerivativeSpace{1},dP2)
        @test_throws MethodError convert(BSplineDerivativeSpace{2},dP2)
    end

    @testset "degenerate" begin
        P1 = BSplineSpace{2}(KnotVector(1:8))
        P2 = BSplineSpace{2}(knotvector"111111151111")
        dP1 = BSplineDerivativeSpace{1}(P1)
        dP2 = BSplineDerivativeSpace{1}(P2)
        @test isdegenerate_R(P1) == isdegenerate_R(dP1)
        @test isdegenerate_R(P2) == isdegenerate_R(dP2)
        @test isdegenerate_I(P1) == isdegenerate_I(dP1)
        @test isdegenerate_I(P2) == isdegenerate_I(dP2)
        for i in 1:dim(P1)
            @test isdegenerate_R(P1, i) == isdegenerate_R(dP1, i)
            @test isdegenerate_I(P1, i) == isdegenerate_I(dP1, i)
        end
        for i in 1:dim(P2)
            @test isdegenerate_R(P2, i) == isdegenerate_R(dP2, i)
            @test isdegenerate_I(P2, i) == isdegenerate_I(dP2, i)
        end
    end

    # Not sure why this @testset doesn't work fine.
    # @testset "$(p)-th degree basis" for p in 0:p_max
    for p in 0:p_max
        k = KnotVector(rand(20)) + (p+1)*KnotVector([0,1])
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
            @test exactdim_R(Pa) == exactdim_R(P)-r
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
        k = KnotVector(rand(10).-1) + KnotVector(rand(10)) + KnotVector(rand(10).+1)
        ts = rand(10)
        for p in 0:5
            P = BSplineSpace{p}(k)
            Q = P
            for r in 0:5
                dP = BSplineDerivativeSpace{r}(P)
                @test (dP == Q) || (r == 0)
                Q = BasicBSpline.derivative(Q)
                for t in ts
                    j = intervalindex(dP,t)
                    B = collect(bsplinebasisall(dP,j,t))
                    @test bsplinebasis.(dP,j:j+p,t) ≈ B
                    @test bsplinebasis₊₀.(dP,j:j+p,t) ≈ B
                    @test bsplinebasis₋₀.(dP,j:j+p,t) ≈ B
                end
            end
        end
    end

    @testset "bsplinebasisall (uniform)" begin
        k = UniformKnotVector(1:20)
        for p in 0:5
            P = BSplineSpace{p}(k)
            for r in 0:5
                dP = BSplineDerivativeSpace{r}(P)
                for _ in 1:10
                    t = rand(domain(dP))
                    j = intervalindex(dP,t)
                    B = collect(bsplinebasisall(dP,j,t))
                    @test bsplinebasis.(dP,j:j+p,t) ≈ B
                    @test bsplinebasis₊₀.(dP,j:j+p,t) ≈ B
                    @test bsplinebasis₋₀.(dP,j:j+p,t) ≈ B
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

        @test bsplinebasis(dP,1,11//5) ===
        bsplinebasis₊₀(dP,1,11//5) ===
        bsplinebasis₋₀(dP,1,11//5) ===
        bsplinebasis′(P,1,11//5)   ===
        bsplinebasis′₊₀(P,1,11//5) ===
        bsplinebasis′₋₀(P,1,11//5) === 16//25
    end

    @testset "Check type" begin
        k = KnotVector{Int}(1:12)
        P0 = BSplineSpace{0}(k)
        P1 = BSplineSpace{1}(k)
        P2 = BSplineSpace{2}(k)

        for r in 0:2
            dP0 = BSplineDerivativeSpace{r}(P0)
            dP1 = BSplineDerivativeSpace{r}(P1)
            dP2 = BSplineDerivativeSpace{r}(P2)
            @test bsplinebasis(dP0,1,5) isa Float64
            @test bsplinebasis(dP1,1,5) isa Float64
            @test bsplinebasis(dP2,1,5) isa Float64

            @test bsplinebasis₊₀(dP0,1,5) isa Float64
            @test bsplinebasis₊₀(dP1,1,5) isa Float64
            @test bsplinebasis₊₀(dP2,1,5) isa Float64

            @test bsplinebasis₋₀(dP0,1,5) isa Float64
            @test bsplinebasis₋₀(dP1,1,5) isa Float64
            @test bsplinebasis₋₀(dP2,1,5) isa Float64

            @test bsplinebasisall(dP0,1,5) isa SVector{1,Float64}
            @test bsplinebasisall(dP1,1,5) isa SVector{2,Float64}
            @test bsplinebasisall(dP2,1,5) isa SVector{3,Float64}
        end
    end

    @testset "Endpoints" begin
        p = 2

        k = KnotVector(1:2) + (p+1)*KnotVector([0,3])
        P0 = BSplineSpace{0}(k)
        P1 = BSplineSpace{1}(k)
        P2 = BSplineSpace{2}(k)

        @test isdegenerate(P0)
        @test isdegenerate(P1)
        @test isnondegenerate(P2)

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
