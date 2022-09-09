function ≒(a,b)
    if a ≈ b
        return true
    elseif abs(a-b) < 1e-10
        return true
    else
        return false
    end
end

@testset "BSplineBasis" begin

    Δt = 1.0e-8
    ε = 1.0e-8
    ∞ = Inf

    @testset "$(p)-th degree basis" for p in 0:4
        Random.seed!(42)

        v = rand(10)
        k = KnotVector(v) + (p+1)*KnotVector([0, 1]) + KnotVector(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == p

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for t in k
            s = sum([bsplinebasis(P, i, t) for i in 1:n])
            s₊₀ = sum([bsplinebasis₊₀(P, i, t) for i in 1:n])
            s₋₀ = sum([bsplinebasis₋₀(P, i, t) for i in 1:n])
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
    end

    @testset "bsplinebasisall" begin
        Random.seed!(42)

        k = KnotVector(rand(10).-1) + KnotVector(rand(10)) + KnotVector(rand(10).+1)
        ts = rand(10)

        for p in 0:5
            P = BSplineSpace{p}(k)
            for t in ts
                j = intervalindex(P,t)
                B = collect(bsplinebasisall(P,j,t))

                _B = bsplinebasis.(P,j:j+p,t)
                @test _B ≈ B

                _B = bsplinebasis₊₀.(P,j:j+p,t)
                @test _B ≈ B

                _B = bsplinebasis₋₀.(P,j:j+p,t)
                @test _B ≈ B
            end
        end
    end

    @testset "Rational" begin
        k = KnotVector{Int}(1:12)
        P = BSplineSpace{3}(k)
        @test k isa KnotVector{Int}
        @test P isa BSplineSpace{3,Int}
        bsplinebasis(P,1,11//5) isa Rational{Int}
        bsplinebasis₊₀(P,1,11//5) isa Rational{Int}
        bsplinebasis₋₀(P,1,11//5) isa Rational{Int}

        @test bsplinebasis(P,1,11//5) ===
        bsplinebasis₊₀(P,1,11//5) ===
        bsplinebasis₋₀(P,1,11//5) === 106//375
    end

    @testset "Check type" begin
        k = KnotVector{Int}(1:12)
        P0 = BSplineSpace{0}(k)
        P1 = BSplineSpace{1}(k)
        P2 = BSplineSpace{2}(k)

        @test bsplinebasis(P0,1,5) isa Float64
        @test bsplinebasis(P1,1,5) isa Float64
        @test bsplinebasis(P2,1,5) isa Float64

        @test bsplinebasis₊₀(P0,1,5) isa Float64
        @test bsplinebasis₊₀(P1,1,5) isa Float64
        @test bsplinebasis₊₀(P2,1,5) isa Float64

        @test bsplinebasis₋₀(P0,1,5) isa Float64
        @test bsplinebasis₋₀(P1,1,5) isa Float64
        @test bsplinebasis₋₀(P2,1,5) isa Float64

        @test bsplinebasisall(P0,1,5) isa SVector{1,Float64}
        @test bsplinebasisall(P1,1,5) isa SVector{2,Float64}
        @test bsplinebasisall(P2,1,5) isa SVector{3,Float64}
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

        @test bsplinebasis.(P0,1:n0,0) == bsplinebasis₊₀.(P0,1:n0,0) == [0,0,1,0,0,0,0]
        @test bsplinebasis.(P1,1:n1,0) == bsplinebasis₊₀.(P1,1:n1,0) == [0,1,0,0,0,0]
        @test bsplinebasis.(P2,1:n2,0) == bsplinebasis₊₀.(P2,1:n2,0) == [1,0,0,0,0]
        @test bsplinebasis.(P0,1:n0,3) == bsplinebasis₋₀.(P0,1:n0,3) == [0,0,0,0,1,0,0]
        @test bsplinebasis.(P1,1:n1,3) == bsplinebasis₋₀.(P1,1:n1,3) == [0,0,0,0,1,0]
        @test bsplinebasis.(P2,1:n2,3) == bsplinebasis₋₀.(P2,1:n2,3) == [0,0,0,0,1]

        @test bsplinebasis₋₀.(P0,1:n0,0) == [0,0,0,0,0,0,0]
        @test bsplinebasis₋₀.(P1,1:n1,0) == [0,0,0,0,0,0]
        @test bsplinebasis₋₀.(P2,1:n2,0) == [0,0,0,0,0]
        @test bsplinebasis₊₀.(P0,1:n0,3) == [0,0,0,0,0,0,0]
        @test bsplinebasis₊₀.(P1,1:n1,3) == [0,0,0,0,0,0]
        @test bsplinebasis₊₀.(P2,1:n2,3) == [0,0,0,0,0]
    end
end
