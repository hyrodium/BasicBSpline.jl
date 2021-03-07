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
    ε = 1.0e-7
    ∞ = Inf

    @testset "0th degree basis" begin
        Random.seed!(42)

        p = 0
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 0

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "1st degree basis" begin
        Random.seed!(42)

        p = 1
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 1

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "2nd degree basis" begin
        Random.seed!(42)

        p = 2
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 2

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "3rd degree basis" begin
        Random.seed!(42)

        p = 3
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 3

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "4th degree basis" begin
        Random.seed!(42)

        p = 4
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 4

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "5th degree basis" begin
        Random.seed!(42)

        p = 5
        v = rand(10)
        k = Knots(v) + (p+1)*Knots(0, 1) + Knots(v[2], v[3])
        P = BSplineSpace(p, k)
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 5

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(i, P, t)
                     ≈ bsplinebasis₊₀(i, fP, t)
                     ≈ bsplinebasis₋₀(i, P, t)
                     ≈ bsplinebasis₋₀(i, fP, t)
                     ≈ bsplinebasis(i, P, t)
                     ≈ bsplinebasis(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(i, P, t) ≈ bsplinebasis₊₀(i, fP, t) ≒ bsplinebasis₊₀(i, P, t₊₀) ≒ bsplinebasis₊₀(i, fP, t₊₀)
                @test bsplinebasis₋₀(i, P, t) ≈ bsplinebasis₋₀(i, fP, t) ≒ bsplinebasis₋₀(i, P, t₋₀) ≒ bsplinebasis₋₀(i, fP, t₋₀)
                @test bsplinebasis(i, P, t) ≈ bsplinebasis(i, fP, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(i, P, t)
                     ≈ bsplinebasis′₊₀(i, fP, t)
                     ≈ bsplinebasis′₋₀(i, P, t)
                     ≈ bsplinebasis′₋₀(i, fP, t)
                     ≈ bsplinebasis′(i, P, t)
                     ≈ bsplinebasis′(i, fP, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(i, P, t) ≈ bsplinebasis′₊₀(i, fP, t) ≒ bsplinebasis′₊₀(i, P, t₊₀) ≒ bsplinebasis′₊₀(i, fP, t₊₀)
                @test bsplinebasis′₋₀(i, P, t) ≈ bsplinebasis′₋₀(i, fP, t) ≒ bsplinebasis′₋₀(i, P, t₋₀) ≒ bsplinebasis′₋₀(i, fP, t₋₀)
                @test bsplinebasis′(i, P, t) ≈ bsplinebasis′(i, fP, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(i, fP, t)
                b′₋₀ = bsplinebasis′₋₀(i, fP, t)
                b′   = bsplinebasis′(i, fP, t)
                a′₊₀ = (bsplinebasis₊₀(i, fP, t+Δt) - bsplinebasis₊₀(i, fP, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(i, fP, t+Δt) - bsplinebasis₋₀(i, fP, t-Δt)) / 2Δt
                a′   = (bsplinebasis(i, fP, t+Δt) - bsplinebasis(i, fP, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
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
