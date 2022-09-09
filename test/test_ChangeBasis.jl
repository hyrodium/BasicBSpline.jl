@testset "ChangeBasis" begin
    ε = 1e-14
    Random.seed!(42)

    function test_changebasis_R(P,P′)
        @test P ⊆ P′
        A = @inferred changebasis(P,P′)
        @test A == BasicBSpline._changebasis_R(P,P′)
        n = dim(P)
        n′ = dim(P′)
        @test size(A) == (n,n′)
        t_min = minimum(knotvector(P)+knotvector(P′))
        t_max = maximum(knotvector(P)+knotvector(P′))
        ts = range(t_min,t_max,length=20)
        for t in ts
            @test norm(bsplinebasis.(P,1:n,t) - A*bsplinebasis.(P′,1:n′,t), Inf) < 1e-14
        end
    end

    @testset "subseteq (R)" begin
        P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
        P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
        P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
        test_changebasis_R(P1, P2)
        test_changebasis_R(P1, P3)
        @test P2 ⊈ P3
    end

    @testset "sqsubseteq (I)" begin
        p4 = 1
        p5 = 2
        P4 = BSplineSpace{p4}(KnotVector([1, 2, 3, 4, 5]))
        P5 = BSplineSpace{p5}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
        @test P4 ⊑ P4
        @test P4 ⊑ P5
        @test P5 ⊒ P4
        @test P5 ⋢ P4
        @test P4 ⋣ P5

        n4, n5 = dim(P4), dim(P5)
        A45 = @inferred changebasis(P4, P5)
        @test size(A45) == (n4, n5)
        Δ45 = [bsplinebasis(P4, i, t) - sum(A45[i, j] * bsplinebasis(P5, j, t) for j in 1:n5) for i in 1:n4, t in 2 * rand(10) .+ 2]
        @test norm(Δ45) < ε

        P6 = BSplineSpace{p4-1}(knotvector(P4)[2:end-1])
        P7 = BSplineSpace{p5-1}(knotvector(P5)[2:end-1])
        @test P6 ⊑ P7
        @test P6 ⊈ P7

        n6, n7 = dim(P6), dim(P7)
        A67 = @inferred changebasis(P6, P7)
        @test size(A67) == (n6, n7)
        Δ67 = [bsplinebasis(P6, i, t) - sum(A67[i, j] * bsplinebasis(P7, j, t) for j in 1:n7) for i in 1:n6, t in 2 * rand(10) .+ 2]
        @test norm(Δ67) < ε
    end

    @testset "changebasis_sim" begin
        for p in 1:3, L in 1:8
            k1 = UniformKnotVector(1:L+2p+1)
            k2 = UniformKnotVector(p+1:p+L+1) + p*KnotVector([p+1]) + KnotVector(p+L+1:2p+L)
            P1 = UniformBSplineSpace{p}(k1)
            P2 = BSplineSpace{p}(k2)
            @test domain(P1) == domain(P2)
            @test dim(P1) == dim(P2)
            @test P1 ≃ P2
            A = @inferred BasicBSpline._changebasis_sim(P1,P2)
            D = domain(P1)
            n = dim(P1)
            for _ in 1:100
                t = rand(D)
                @test norm(bsplinebasis.(P1,1:n,t) - A*bsplinebasis.(P2,1:n,t), Inf) < 1e-14
            end
        end
    end

    @testset "non-float" begin
        k = UniformKnotVector(1:15)
        P = UniformBSplineSpace{3}(k)
        @test changebasis(P,P) isa Matrix{Float64}
        @test BasicBSpline._changebasis_R(P,P) isa Matrix{Float64}
        @test BasicBSpline._changebasis_I(P,P) isa Matrix{Float64}
        @test BasicBSpline._changebasis_sim(P,P) isa Matrix{Float64}

        k = UniformKnotVector(1:15//1)
        P = UniformBSplineSpace{3}(k)
        @test changebasis(P,P) isa Matrix{Rational{Int}}
        @test BasicBSpline._changebasis_R(P,P) isa Matrix{Rational{Int}}
        @test BasicBSpline._changebasis_I(P,P) isa Matrix{Rational{Int}}
        @test BasicBSpline._changebasis_sim(P,P) isa Matrix{Rational{Int}}

        k = UniformKnotVector(1:BigInt(15))
        P = UniformBSplineSpace{3}(k)
        @test changebasis(P,P) isa Matrix{BigFloat}
        @test BasicBSpline._changebasis_R(P,P) isa Matrix{BigFloat}
        @test BasicBSpline._changebasis_I(P,P) isa Matrix{BigFloat}
        @test BasicBSpline._changebasis_sim(P,P) isa Matrix{BigFloat}
    end

    @testset "uniform" begin
        for r in 1:5
            k = UniformKnotVector(0:r:50)
            for p in 0:4
                P = UniformBSplineSpace{p}(k)

                k′ = UniformKnotVector(0:50)
                P′ = UniformBSplineSpace{p}(k′)
                A1 = @inferred changebasis(P, P′)
                A2 = @inferred changebasis(BSplineSpace(P), BSplineSpace(P′))
                @test P ⊆ P′
                @test A1 ≈ A2

                left = leftendpoint(domain(P))-p
                right = rightendpoint(domain(P))+p
                k′ = UniformKnotVector(left:right)
                P′ = UniformBSplineSpace{p}(k′)
                A3 = @inferred changebasis(P, P′)
                A4 = @inferred changebasis(BSplineSpace(P), BSplineSpace(P′))
                @test P ⊑ P′
                @test A3 ≈ A4
            end
        end
    end

    @testset "derivative" begin
        k1 = KnotVector(5*rand(12))
        k2 = k1 + KnotVector(rand(3)*3 .+ 1)
        p = 5
        dP1 = BSplineDerivativeSpace{1}(BSplineSpace{p}(k1))
        dP2 = BSplineDerivativeSpace{0}(BSplineSpace{p-1}(k1))
        P1 = BSplineSpace{p-1}(k1)
        P2 = BSplineSpace{p-1}(k2)

        test_changebasis_R(dP1, P1)
        test_changebasis_R(dP1, P2)
        test_changebasis_R(dP2, P1)
        test_changebasis_R(dP2, P2)
        test_changebasis_R(dP2, P1)
        test_changebasis_R(dP2, P2)

        test_changebasis_R(dP1, dP2)
        test_changebasis_R(dP1, P2)
        test_changebasis_R(dP2, dP2)
        test_changebasis_R(dP2, P2)
        test_changebasis_R(dP2, dP2)
        test_changebasis_R(dP2, P2)

        dP2 = BSplineDerivativeSpace{2}(BSplineSpace{p}(k1))
        dP3 = BSplineDerivativeSpace{1}(BSplineSpace{p-1}(k1))
        dP4 = BSplineDerivativeSpace{0}(BSplineSpace{p-2}(k1))
        dP5 = BSplineDerivativeSpace{2}(BSplineSpace{p}(k2))
        P3 = BSplineSpace{p-2}(k1)
        P4 = BSplineSpace{p-2}(k2)

        test_changebasis_R(dP2, dP3)
        @test dP2 ⊉ dP3
        test_changebasis_R(dP3, dP4)
        @test dP3 ⊉ dP4
        test_changebasis_R(dP4, P3)
        test_changebasis_R(P3, dP4)

        test_changebasis_R(dP2, dP5)
        @test dP2 ⊉ dP5
        test_changebasis_R(dP5, dP5)
        test_changebasis_R(dP2, dP2)
        test_changebasis_R(dP5, P4)
        @test dP5 ⊉ P4
        test_changebasis_R(dP3, P4)
        @test dP3 ⊉ P4

        @test_throws DomainError changebasis(P4, P3)
    end
end
