@testset "RationalBSplineManifold" begin
    @testset "1dim-arc" begin
        a = [SVector(1,0), SVector(1,1), SVector(0,1), SVector(-1,1), SVector(-1,0)]
        w = [1, 1/√2, 1, 1/√2, 1]

        k = KnotVector([0,0,0,1,1,2,2,2])
        p = 2
        P = BSplineSpace{p}(k)

        M = RationalBSplineManifold(a,w,(P,))
        for _ in 1:100
            t1 = 2rand()
            @test norm(M(t1)) ≈ 1
        end
        @test_throws DomainError M(-0.1)
        @test_throws DomainError M(2.1)
    end

    @testset "compatibility with BSplineManifold" begin
        k1 = KnotVector(rand(30)) + 5KnotVector(0,1)
        k2 = KnotVector(rand(30)) + 5KnotVector(0,1)
        k3 = KnotVector(rand(30)) + 5KnotVector(0,1)
        p1 = 2
        p2 = 3
        p3 = 4
        P1 = BSplineSpace{p1}(k1)
        P2 = BSplineSpace{p2}(k2)
        P3 = BSplineSpace{p3}(k3)
        n1 = dim(P1)
        n2 = dim(P2)
        n3 = dim(P3)

        @testset "1dim" begin
            a = randn(n1)
            w1 = ones(n1)
            w2 = 2ones(n1)
            M = BSplineManifold(a,(P1,))
            M1 = RationalBSplineManifold(a,w1,(P1,))
            M2 = RationalBSplineManifold(a,w2,(P1,))
            for _ in 1:100
                t1 = rand()
                @test M(t1) ≈ M1(t1) atol=1e-14
                @test M(t1) ≈ M2(t1) atol=1e-14
            end
        end

        @testset "2dim" begin
            a = randn(n1,n2)
            w1 = ones(n1,n2)
            w2 = 2ones(n1,n2)
            M = BSplineManifold(a,(P1,P2))
            M1 = RationalBSplineManifold(a,w1,(P1,P2))
            M2 = RationalBSplineManifold(a,w2,(P1,P2))
            for _ in 1:100
                t1 = rand()
                t2 = rand()
                @test M(t1,t2) ≈ M1(t1,t2) atol=1e-14
                @test M(t1,t2) ≈ M2(t1,t2) atol=1e-14
            end
        end

        @testset "3dim" begin
            a = randn(n1,n2,n3)
            w1 = ones(n1,n2,n3)
            w2 = 2ones(n1,n2,n3)
            M = BSplineManifold(a,(P1,P2,P3))
            M1 = RationalBSplineManifold(a,w1,(P1,P2,P3))
            M2 = RationalBSplineManifold(a,w2,(P1,P2,P3))
            for _ in 1:100
                t1 = rand()
                t2 = rand()
                t3 = rand()
                @test M(t1,t2,t3) ≈ M1(t1,t2,t3) atol=1e-14
                @test M(t1,t2,t3) ≈ M2(t1,t2,t3) atol=1e-14
            end
        end
    end
end
