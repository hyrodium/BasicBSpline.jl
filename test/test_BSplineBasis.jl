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
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 0

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₊₀(fP, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis₋₀(fP, i, t)
                     ≈ bsplinebasis(P, i, t)
                     ≈ bsplinebasis(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≈ bsplinebasis₊₀(fP, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀) ≒ bsplinebasis₊₀(fP, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≈ bsplinebasis₋₀(fP, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀) ≒ bsplinebasis₋₀(fP, i, t₋₀)
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(fP, i, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₊₀(fP, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′₋₀(fP, i, t)
                     ≈ bsplinebasis′(P, i, t)
                     ≈ bsplinebasis′(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≈ bsplinebasis′₊₀(fP, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀) ≒ bsplinebasis′₊₀(fP, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≈ bsplinebasis′₋₀(fP, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀) ≒ bsplinebasis′₋₀(fP, i, t₋₀)
                @test bsplinebasis′(P, i, t) ≈ bsplinebasis′(fP, i, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(fP, i, t)
                b′₋₀ = bsplinebasis′₋₀(fP, i, t)
                b′   = bsplinebasis′(fP, i, t)
                a′₊₀ = (bsplinebasis₊₀(fP, i, t+Δt) - bsplinebasis₊₀(fP, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(fP, i, t+Δt) - bsplinebasis₋₀(fP, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(fP, i, t+Δt) - bsplinebasis(fP, i, t-Δt)) / 2Δt
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
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 1

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₊₀(fP, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis₋₀(fP, i, t)
                     ≈ bsplinebasis(P, i, t)
                     ≈ bsplinebasis(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≈ bsplinebasis₊₀(fP, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀) ≒ bsplinebasis₊₀(fP, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≈ bsplinebasis₋₀(fP, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀) ≒ bsplinebasis₋₀(fP, i, t₋₀)
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(fP, i, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₊₀(fP, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′₋₀(fP, i, t)
                     ≈ bsplinebasis′(P, i, t)
                     ≈ bsplinebasis′(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≈ bsplinebasis′₊₀(fP, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀) ≒ bsplinebasis′₊₀(fP, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≈ bsplinebasis′₋₀(fP, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀) ≒ bsplinebasis′₋₀(fP, i, t₋₀)
                @test bsplinebasis′(P, i, t) ≈ bsplinebasis′(fP, i, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(fP, i, t)
                b′₋₀ = bsplinebasis′₋₀(fP, i, t)
                b′   = bsplinebasis′(fP, i, t)
                a′₊₀ = (bsplinebasis₊₀(fP, i, t+Δt) - bsplinebasis₊₀(fP, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(fP, i, t+Δt) - bsplinebasis₋₀(fP, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(fP, i, t+Δt) - bsplinebasis(fP, i, t-Δt)) / 2Δt
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
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 2

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₊₀(fP, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis₋₀(fP, i, t)
                     ≈ bsplinebasis(P, i, t)
                     ≈ bsplinebasis(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≈ bsplinebasis₊₀(fP, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀) ≒ bsplinebasis₊₀(fP, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≈ bsplinebasis₋₀(fP, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀) ≒ bsplinebasis₋₀(fP, i, t₋₀)
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(fP, i, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₊₀(fP, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′₋₀(fP, i, t)
                     ≈ bsplinebasis′(P, i, t)
                     ≈ bsplinebasis′(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≈ bsplinebasis′₊₀(fP, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀) ≒ bsplinebasis′₊₀(fP, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≈ bsplinebasis′₋₀(fP, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀) ≒ bsplinebasis′₋₀(fP, i, t₋₀)
                @test bsplinebasis′(P, i, t) ≈ bsplinebasis′(fP, i, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(fP, i, t)
                b′₋₀ = bsplinebasis′₋₀(fP, i, t)
                b′   = bsplinebasis′(fP, i, t)
                a′₊₀ = (bsplinebasis₊₀(fP, i, t+Δt) - bsplinebasis₊₀(fP, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(fP, i, t+Δt) - bsplinebasis₋₀(fP, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(fP, i, t+Δt) - bsplinebasis(fP, i, t-Δt)) / 2Δt
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
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 3

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₊₀(fP, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis₋₀(fP, i, t)
                     ≈ bsplinebasis(P, i, t)
                     ≈ bsplinebasis(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≈ bsplinebasis₊₀(fP, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀) ≒ bsplinebasis₊₀(fP, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≈ bsplinebasis₋₀(fP, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀) ≒ bsplinebasis₋₀(fP, i, t₋₀)
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(fP, i, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₊₀(fP, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′₋₀(fP, i, t)
                     ≈ bsplinebasis′(P, i, t)
                     ≈ bsplinebasis′(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≈ bsplinebasis′₊₀(fP, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀) ≒ bsplinebasis′₊₀(fP, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≈ bsplinebasis′₋₀(fP, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀) ≒ bsplinebasis′₋₀(fP, i, t₋₀)
                @test bsplinebasis′(P, i, t) ≈ bsplinebasis′(fP, i, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(fP, i, t)
                b′₋₀ = bsplinebasis′₋₀(fP, i, t)
                b′   = bsplinebasis′(fP, i, t)
                a′₊₀ = (bsplinebasis₊₀(fP, i, t+Δt) - bsplinebasis₊₀(fP, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(fP, i, t+Δt) - bsplinebasis₋₀(fP, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(fP, i, t+Δt) - bsplinebasis(fP, i, t-Δt)) / 2Δt
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
        fP = FastBSplineSpace(p, k)
        n = dim(P)

        @test degree(P) == degree(fP) == 4

        # Check the values of B-spline basis funcitons
        for i in 1:n, t in rand(5)
            for t in rand(5)
                @test (bsplinebasis₊₀(P, i, t)
                     ≈ bsplinebasis₊₀(fP, i, t)
                     ≈ bsplinebasis₋₀(P, i, t)
                     ≈ bsplinebasis₋₀(fP, i, t)
                     ≈ bsplinebasis(P, i, t)
                     ≈ bsplinebasis(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis₊₀(P, i, t) ≈ bsplinebasis₊₀(fP, i, t) ≒ bsplinebasis₊₀(P, i, t₊₀) ≒ bsplinebasis₊₀(fP, i, t₊₀)
                @test bsplinebasis₋₀(P, i, t) ≈ bsplinebasis₋₀(fP, i, t) ≒ bsplinebasis₋₀(P, i, t₋₀) ≒ bsplinebasis₋₀(fP, i, t₋₀)
                @test bsplinebasis(P, i, t) ≈ bsplinebasis(fP, i, t)
            end
        end

        # Check the values of the derivative of B-spline basis funcitons
        for i in 1:n
            for t in rand(5)
                @test (bsplinebasis′₊₀(P, i, t)
                     ≈ bsplinebasis′₊₀(fP, i, t)
                     ≈ bsplinebasis′₋₀(P, i, t)
                     ≈ bsplinebasis′₋₀(fP, i, t)
                     ≈ bsplinebasis′(P, i, t)
                     ≈ bsplinebasis′(fP, i, t))
            end
            for t in k
                t₊₀ = nextfloat(t)
                t₋₀ = prevfloat(t)
                @test bsplinebasis′₊₀(P, i, t) ≈ bsplinebasis′₊₀(fP, i, t) ≒ bsplinebasis′₊₀(P, i, t₊₀) ≒ bsplinebasis′₊₀(fP, i, t₊₀)
                @test bsplinebasis′₋₀(P, i, t) ≈ bsplinebasis′₋₀(fP, i, t) ≒ bsplinebasis′₋₀(P, i, t₋₀) ≒ bsplinebasis′₋₀(fP, i, t₋₀)
                @test bsplinebasis′(P, i, t) ≈ bsplinebasis′(fP, i, t)
            end
        end

        # Check the derivative
        @testset "derivative" begin
            Random.seed!(57)
            for i in 1:n, t in rand(5)
                b′₊₀ = bsplinebasis′₊₀(fP, i, t)
                b′₋₀ = bsplinebasis′₋₀(fP, i, t)
                b′   = bsplinebasis′(fP, i, t)
                a′₊₀ = (bsplinebasis₊₀(fP, i, t+Δt) - bsplinebasis₊₀(fP, i, t-Δt)) / 2Δt
                a′₋₀ = (bsplinebasis₋₀(fP, i, t+Δt) - bsplinebasis₋₀(fP, i, t-Δt)) / 2Δt
                a′   = (bsplinebasis(fP, i, t+Δt) - bsplinebasis(fP, i, t-Δt)) / 2Δt
                @test b′₊₀ ≈ b′₋₀ ≈ b′ ≈ a′₊₀ ≈ a′₋₀ ≈ a′
            end
        end
    end

    @testset "_bsplinebasis" begin
        Random.seed!(42)

        k = Knots(rand(10).-1) + Knots(rand(10)) + Knots(rand(10).+1)
        ts = rand(10)

        for p in 0:BasicBSpline.MAX_DEGREE
            fP = FastBSplineSpace(p,k)
            P = BSplineSpace{p}(k)
            for t in ts
                j = BasicBSpline._knotindex(fP,t)
                J = intervalindex(P,t)
                _B = BasicBSpline._bsplinebasis(fP,t,j)
                _b = bsplinebasisall(P,J,t)
                # _B′ = BasicBSpline._bsplinebasis′(fP,t,j)
                B = [bsplinebasis(fP,i,t) for i in j-p:j]
                B′ = [bsplinebasis′(fP,i,t) for i in j-p:j]
                @test norm(collect(_B) - B) < ε

                _B = [bsplinebasis(P,i,t) for i in j-p:j]
                @test norm(_B - B) < ε

                _B = [bsplinebasis₊₀(P,i,t) for i in j-p:j]
                @test norm(_B - B) < ε

                _B = [bsplinebasis₋₀(P,i,t) for i in j-p:j]
                @test norm(_B - B) < ε

                _B′ =[ bsplinebasis′(P,i,t) for i in j-p:j]
                @test norm(_B′ - B′) < ε

                _B′ =[ bsplinebasis′₊₀(P,i,t) for i in j-p:j]
                @test norm(_B′ - B′) < ε

                _B′ =[ bsplinebasis′₋₀(P,i,t) for i in j-p:j]
                @test norm(_B′ - B′) < ε
            end
        end
    end
end
