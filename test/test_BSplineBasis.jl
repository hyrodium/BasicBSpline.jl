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
        k = KnotVector(v) + (p+1)*KnotVector(0, 1) + KnotVector(v[2], v[3])
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
                _B = [bsplinebasis(P,i,t) for i in j:j+p]
                @test _B ≈ B

                _B = [bsplinebasis(P,i,t) for i in j:j+p]
                @test _B ≈ B

                _B = [bsplinebasis₊₀(P,i,t) for i in j:j+p]
                @test _B ≈ B

                _B = [bsplinebasis₋₀(P,i,t) for i in j:j+p]
                @test _B ≈ B
            end
        end
    end

    @testset "Rational" begin
        k = KnotVector{Int}(1:12)
        P = BSplineSpace{3}(k)
        @test k isa KnotVector{Int}
        @test P isa BSplineSpace{3,Int}
        @test bsplinebasis(P,1,11//5) isa Rational{Int}
        @test bsplinebasis(P,1,11//5) === bsplinebasis(P,2,16//5) === 106//375
    end
end
