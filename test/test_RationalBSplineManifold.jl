@testset "RationalBSplineManifold" begin
    @testset "constructor" begin
        k1 = KnotVector(Float64[0,0,1,1])
        k2 = UniformKnotVector(1.0:8.0)
        k3 = KnotVector(rand(8))
        _k1 = KnotVector([0.0, 0.0, 1.0, 1.0])
        _k2 = KnotVector(1.0:8.0)
        _k3 = view(k3,:)
        P1 = BSplineSpace{1}(k1)
        P2 = BSplineSpace{3}(k2)
        P3 = BSplineSpace{2}(k3)
        _P1 = BSplineSpace{1}(_k1)
        _P2 = BSplineSpace{3}(_k2)
        _P3 = BSplineSpace{2}(_k3)

        # 0-dim
        a = fill(1.2)
        w = fill(4.2)
        M1 = RationalBSplineManifold{0,(),Float64,Float64,Int,Tuple{}}(a, w, ())
        M2 = RationalBSplineManifold{0,(),Float64,Float64,Int,Tuple{}}(copy(a), copy(w), ())
        M3 = BSplineManifold{0,(),Float64,Int,Tuple{}}(a, ())
        M4 = RationalBSplineManifold(M3)
        @test M1 == M1
        @test M2 == M2
        @test M3 == M3
        @test M4 == M4
        @test M1 == M2 == RationalBSplineManifold(M1) == RationalBSplineManifold(M2)
        @test hash(M1) == hash(M2) == hash(RationalBSplineManifold(M1)) == hash(RationalBSplineManifold(M2))
        @test M3 ≠ M4 == RationalBSplineManifold(M4)
        @test hash(M3) ≠ hash(M4) == hash(RationalBSplineManifold(M4))
        @test M1() == 1.2
        @test M2() == 1.2
        @test M3() == 1.2
        @test M4() == 1.2
        @test all(isone.(weights(M4)))

        # 4-dim
        a = rand(dim(P1), dim(P2), dim(P3), dim(P3))
        w = rand(dim(P1), dim(P2), dim(P3), dim(P3))
        @test RationalBSplineManifold(a, w, P1, P2, P3, P3) == RationalBSplineManifold(a, w, P1, P2, P3, P3)
        @test hash(RationalBSplineManifold(a, w, P1, P2, P3, P3)) == hash(RationalBSplineManifold(a, w, P1, P2, P3, P3))
        @test RationalBSplineManifold(a, w, P1, P2, P3, P3) == RationalBSplineManifold(a, w, _P1, _P2, _P3, _P3)
        @test hash(RationalBSplineManifold(a, w, P1, P2, P3, P3)) == hash(RationalBSplineManifold(a, w, _P1, _P2, _P3, _P3))
        @test RationalBSplineManifold(a, w, P1, P2, P3, P3) == RationalBSplineManifold(copy(a), copy(w), _P1, _P2, _P3, _P3)
        @test hash(RationalBSplineManifold(a, w, P1, P2, P3, P3)) == hash(RationalBSplineManifold(copy(a), copy(w), _P1, _P2, _P3, _P3))
    end

    @testset "1dim-arc" begin
        a = [SVector(1,0), SVector(1,1), SVector(0,1), SVector(-1,1), SVector(-1,0)]
        w = [1, 1/√2, 1, 1/√2, 1]

        k = KnotVector([0,0,0,1,1,2,2,2])
        p = 2
        P = BSplineSpace{p}(k)

        M = RationalBSplineManifold(a,w,P)
        for _ in 1:100
            t1 = 2rand()
            @test norm(M(t1)) ≈ 1
        end
        @test_throws DomainError M(-0.1)
        @test_throws DomainError M(2.1)
    end

    @testset "general cases" begin
        p1 = 2
        p2 = 3
        p3 = 4
        k1 = KnotVector(rand(30)) + (p1+1)*KnotVector([0,1])
        k2 = KnotVector(rand(30)) + (p2+1)*KnotVector([0,1])
        k3 = KnotVector(rand(30)) + (p3+1)*KnotVector([0,1])
        P1 = BSplineSpace{p1}(k1)
        P2 = BSplineSpace{p2}(k2)
        P3 = BSplineSpace{p3}(k3)
        n1 = dim(P1)
        n2 = dim(P2)
        n3 = dim(P3)

        @testset "RationalBSpilneManifold definition" begin
            @testset "1dim" begin
                a = randn(ComplexF64,n1)
                w = ones(n1)+rand(n1)
                R = RationalBSplineManifold(a,w,(P1,))
                A = BSplineManifold(a.*w,(P1,))
                W = BSplineManifold(w,(P1,))
                @test R(:) == R
                @test Base.mightalias(controlpoints(R), controlpoints(R))
                @test !Base.mightalias(controlpoints(R), controlpoints(R(:)))
                ts = [(rand.(domain.(bsplinespaces(R)))) for _ in 1:10]
                for (t1,) in ts
                    @test R(t1) ≈ A(t1)/W(t1)
                end
            end
            @testset "2dim" begin
                a = randn(ComplexF64,n1,n2)
                w = ones(n1,n2)+rand(n1,n2)
                R = RationalBSplineManifold(a,w,(P1,P2))
                A = BSplineManifold(a.*w,(P1,P2))
                W = BSplineManifold(w,(P1,P2))
                @test R(:,:) == R
                @test Base.mightalias(controlpoints(R), controlpoints(R))
                @test !Base.mightalias(controlpoints(R), controlpoints(R(:,:)))
                ts = [(rand.(domain.(bsplinespaces(R)))) for _ in 1:10]
                for (t1, t2) in ts
                    @test R(t1,t2) ≈ A(t1,t2)/W(t1,t2)
                    @test R(t1,t2) ≈ R(t1,:)(t2)
                    @test R(t1,t2) ≈ R(:,t2)(t1)
                end
            end
            @testset "3dim" begin
                a = randn(ComplexF64,n1,n2,n3)
                w = ones(n1,n2,n3)+rand(n1,n2,n3)
                R = RationalBSplineManifold(a,w,(P1,P2,P3))
                A = BSplineManifold(a.*w,(P1,P2,P3))
                W = BSplineManifold(w,(P1,P2,P3))
                @test R(:,:,:) == R
                @test Base.mightalias(controlpoints(R), controlpoints(R))
                @test !Base.mightalias(controlpoints(R), controlpoints(R(:,:,:)))
                ts = [(rand.(domain.(bsplinespaces(R)))) for _ in 1:10]
                for (t1, t2, t3) in ts
                    @test R(t1,t2,t3) ≈ A(t1,t2,t3)/W(t1,t2,t3)
                    @test R(t1,t2,t3) ≈ R(t1,:,:)(t2,t3)
                    @test R(t1,t2,t3) ≈ R(:,t2,:)(t1,t3)
                    @test R(t1,t2,t3) ≈ R(:,:,t3)(t1,t2)
                    @test R(t1,t2,t3) ≈ R(t1,t2,:)(t3)
                    @test R(t1,t2,t3) ≈ R(t1,:,t3)(t2)
                    @test R(t1,t2,t3) ≈ R(:,t2,t3)(t1)
                end
            end
            @testset "4dim" begin
                @testset "BSplineManifold-4dim" begin
                    P1 = BSplineSpace{3}(knotvector"43211112112112")
                    P2 = BSplineSpace{4}(knotvector"3211  1 112112115")
                    P3 = BSplineSpace{5}(knotvector" 411  1 112112 11133")
                    P4 = BSplineSpace{4}(knotvector"4 1112113 11 13 3")
                    n1 = dim(P1)
                    n2 = dim(P2)
                    n3 = dim(P3)
                    n4 = dim(P4)
                    a = rand(n1, n2, n3, n4)
                    w = rand(n1, n2, n3, n4) .+ 1
                    R = RationalBSplineManifold(a, w, (P1, P2, P3, P4))
                    @test dim(R) == 4

                    ts = [(rand.(domain.((P1, P2, P3, P4)))) for _ in 1:10]
                    for (t1, t2, t3, t4) in ts
                        @test R(t1,t2,t3,t4) ≈ R(t1,:,:,:)(   t2,t3,t4)
                        @test R(t1,t2,t3,t4) ≈ R(:,t2,:,:)(t1,   t3,t4)
                        @test R(t1,t2,t3,t4) ≈ R(:,:,t3,:)(t1,t2,   t4)
                        @test R(t1,t2,t3,t4) ≈ R(:,:,:,t4)(t1,t2,t3   )

                        @test R(t1,t2,t3,t4) ≈ R(t1,t2,:,:)(t3,t4)
                        @test R(t1,t2,t3,t4) ≈ R(t1,:,t3,:)(t2,t4)
                        @test R(t1,t2,t3,t4) ≈ R(t1,:,:,t4)(t2,t3)
                        @test R(t1,t2,t3,t4) ≈ R(:,t2,t3,:)(t1,t4)
                        @test R(t1,t2,t3,t4) ≈ R(:,t2,:,t4)(t1,t3)
                        @test R(t1,t2,t3,t4) ≈ R(:,:,t3,t4)(t1,t2)

                        @test R(t1,t2,t3,t4) ≈ R(:,t2,t3,t4)(t1)
                        @test R(t1,t2,t3,t4) ≈ R(t1,:,t3,t4)(t2)
                        @test R(t1,t2,t3,t4) ≈ R(t1,t2,:,t4)(t3)
                        @test R(t1,t2,t3,t4) ≈ R(t1,t2,t3,:)(t4)
                    end
                    @test_throws DomainError R(-5,-8,-120,1)

                    @test R(:,:,:,:) == R
                    @test Base.mightalias(controlpoints(R), controlpoints(R))
                    @test !Base.mightalias(controlpoints(R), controlpoints(R(:,:,:,:)))
                end
            end
        end

        @testset "compatibility with BSplineManifold" begin
            @testset "1dim" begin
                a = randn(n1)
                w1 = ones(n1)
                w2 = 2ones(n1)
                M = BSplineManifold(a,P1)
                M1 = RationalBSplineManifold(a,w1,P1)
                M2 = RationalBSplineManifold(a,w2,P1)
                ts = [(rand.(domain.(bsplinespaces(M)))) for _ in 1:10]
                for (t1,) in ts
                    @test M(t1) ≈ M1(t1) atol=1e-14
                    @test M(t1) ≈ M2(t1) atol=1e-14
                end
            end
            @testset "2dim" begin
                a = randn(n1,n2)
                w1 = ones(n1,n2)
                w2 = 2ones(n1,n2)
                M = BSplineManifold(a,P1,P2)
                M1 = RationalBSplineManifold(a,w1,P1,P2)
                M2 = RationalBSplineManifold(a,w2,P1,P2)
                ts = [(rand.(domain.(bsplinespaces(M)))) for _ in 1:10]
                for (t1, t2) in ts
                    @test M(t1,t2) ≈ M1(t1,t2) atol=1e-14
                    @test M(t1,t2) ≈ M2(t1,t2) atol=1e-14
                end
            end
            @testset "3dim" begin
                a = randn(n1,n2,n3)
                w1 = ones(n1,n2,n3)
                w2 = 2ones(n1,n2,n3)
                M = BSplineManifold(a,P1,P2,P3)
                M1 = RationalBSplineManifold(a,w1,P1,P2,P3)
                M2 = RationalBSplineManifold(a,w2,P1,P2,P3)
                ts = [(rand.(domain.(bsplinespaces(M)))) for _ in 1:10]
                for (t1, t2, t3) in ts
                    @test M(t1,t2,t3) ≈ M1(t1,t2,t3) atol=1e-14
                    @test M(t1,t2,t3) ≈ M2(t1,t2,t3) atol=1e-14
                end
            end
        end
    end
end
