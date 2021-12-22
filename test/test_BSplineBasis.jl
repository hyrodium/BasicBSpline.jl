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

    @testset "0th degree basis" begin
        Random.seed!(42)

        p = 0
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == 0

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
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(P, i, t)
                b′₋₀ = bsplinebasis′₋₀(P, i, t)
                b′   = bsplinebasis′(P, i, t)
                a′₊₀ = (bsplinebasis₊₀(P, i, t+Δt) - bsplinebasis₊₀(P, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(P, i, t+Δt) - bsplinebasis₋₀(P, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(P, i, t+Δt) - bsplinebasis(P, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "1st degree basis" begin
        Random.seed!(42)

        p = 1
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == 1

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
                @test bsplinebasis₋₀(P, i, t)≒ bsplinebasis₋₀(P, i, t₋₀)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(P, i, t)
                b′₋₀ = bsplinebasis′₋₀(P, i, t)
                b′   = bsplinebasis′(P, i, t)
                a′₊₀ = (bsplinebasis₊₀(P, i, t+Δt) - bsplinebasis₊₀(P, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(P, i, t+Δt) - bsplinebasis₋₀(P, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(P, i, t+Δt) - bsplinebasis(P, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "2nd degree basis" begin
        Random.seed!(42)

        p = 2
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == 2

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
                @test bsplinebasis₋₀(P, i, t)≒ bsplinebasis₋₀(P, i, t₋₀)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(P, i, t)
                b′₋₀ = bsplinebasis′₋₀(P, i, t)
                b′   = bsplinebasis′(P, i, t)
                a′₊₀ = (bsplinebasis₊₀(P, i, t+Δt) - bsplinebasis₊₀(P, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(P, i, t+Δt) - bsplinebasis₋₀(P, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(P, i, t+Δt) - bsplinebasis(P, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "3rd degree basis" begin
        Random.seed!(42)

        p = 3
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == 3

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
                @test bsplinebasis₋₀(P, i, t)≒ bsplinebasis₋₀(P, i, t₋₀)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(P, i, t)
                b′₋₀ = bsplinebasis′₋₀(P, i, t)
                b′   = bsplinebasis′(P, i, t)
                a′₊₀ = (bsplinebasis₊₀(P, i, t+Δt) - bsplinebasis₊₀(P, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(P, i, t+Δt) - bsplinebasis₋₀(P, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(P, i, t+Δt) - bsplinebasis(P, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "4th degree basis" begin
        Random.seed!(42)

        p = 4
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace{p}(k)
        n = dim(P)

        @test degree(P) == 4

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
                @test bsplinebasis₋₀(P, i, t)≒ bsplinebasis₋₀(P, i, t₋₀)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′(P, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(P, i, t)
                b′₋₀ = bsplinebasis′₋₀(P, i, t)
                b′   = bsplinebasis′(P, i, t)
                a′₊₀ = (bsplinebasis₊₀(P, i, t+Δt) - bsplinebasis₊₀(P, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(P, i, t+Δt) - bsplinebasis₋₀(P, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(P, i, t+Δt) - bsplinebasis(P, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "_bsplinebasis" begin
        Random.seed!(42)

        k = Knots(rand(10).-1) + Knots(rand(10)) + Knots(rand(10).+1)
        ts = rand(10)

        for p in 0:BasicBSpline.MAX_DEGREE
            P = BSplineSpace{p}(k)
            for t in ts
                j = intervalindex(P,t)
                _B = bsplinebasisall(P,j,t)
                B = [bsplinebasis(P,i,t) for i in j:j+p]
                B′ = [bsplinebasis′(P,i,t) for i in j:j+p]
                @test norm(collect(_B) - B) < ε

                _B = [bsplinebasis(P,i,t) for i in j:j+p]
                @test norm(_B - B) < ε

                _B = [bsplinebasis₊₀(P,i,t) for i in j:j+p]
                @test norm(_B - B) < ε

                _B = [bsplinebasis₋₀(P,i,t) for i in j:j+p]
                @test norm(_B - B) < ε

                _B′ =[ bsplinebasis′(P,i,t) for i in j:j+p]
                @test norm(_B′ - B′) < ε

                _B′ =[ bsplinebasis′₊₀(P,i,t) for i in j:j+p]
                @test norm(_B′ - B′) < ε

                _B′ =[ bsplinebasis′₋₀(P,i,t) for i in j:j+p]
                @test norm(_B′ - B′) < ε
            end
        end
    end
end
