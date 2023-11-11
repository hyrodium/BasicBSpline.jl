@testset "ChangeBasis" begin
    ε = 1e-13
    function test_changebasis_R(P,P′)
        @test P ⊆ P′
        A = @inferred changebasis(P,P′)
        @test A isa SparseMatrixCSC
        @test !any(iszero.(A.nzval))
        @test A == BasicBSpline._changebasis_R(P,P′) == changebasis_R(P,P′)
        n = dim(P)
        n′ = dim(P′)
        @test size(A) == (n,n′)
        ts = range(extrema(knotvector(P)+knotvector(P′))..., length=20)
        for t in ts
            @test norm(bsplinebasis.(P,1:n,t) - A*bsplinebasis.(P′,1:n′,t), Inf) < ε
        end
        @test iszero(view(A, findall(BasicBSpline._iszeros_R(P)), :))
        @test iszero(view(A, :, findall(BasicBSpline._iszeros_R(P′))))
    end

    function test_changebasis_I(P1,P2)
        @test P1 ⊑ P2
        A = @inferred changebasis_I(P1,P2)
        @test A isa SparseMatrixCSC
        @test !any(iszero.(A.nzval))
        n1 = dim(P1)
        n2 = dim(P2)
        @test size(A) == (n1,n2)
        d = domain(P1)
        ts = range(extrema(d)..., length=21)[2:end-1]
        for t in ts
            @test norm(bsplinebasis.(P1,1:n1,t) - A*bsplinebasis.(P2,1:n2,t), Inf) < ε
        end
    end

    @testset "subseteq (R)" begin
        P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
        P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
        P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
        P4 = BSplineSpace{1}(KnotVector([1, 3, 4, 4, 4, 4, 5, 8]))
        P5 = expandspace_R(P3, Val(2))
        P6 = expandspace_R(P3, KnotVector([1.2]))
        P7 = expandspace_R(P3, KnotVector([1, 1.2]))

        test_changebasis_R(P1, P2)
        test_changebasis_R(P1, P3)
        test_changebasis_R(P1, P4)
        test_changebasis_R(P1, P1)
        test_changebasis_R(P2, P2)
        test_changebasis_R(P3, P3)
        test_changebasis_R(P4, P4)

        test_changebasis_R(P1, P3)
        test_changebasis_R(P3, P5)
        test_changebasis_R(P1, P5)
        A13 = changebasis(P1, P3)
        A35 = changebasis(P3, P5)
        A15 = changebasis(P1, P5)
        @test A15 ≈ A13 * A35

        test_changebasis_R(P1, P3)
        test_changebasis_R(P3, P6)
        test_changebasis_R(P1, P6)
        A13 = changebasis(P1, P3)
        A36 = changebasis(P3, P6)
        A16 = changebasis(P1, P6)
        A16 ≈ A13 * A36

        test_changebasis_R(P3, P6)
        test_changebasis_R(P6, P7)
        test_changebasis_R(P3, P7)
        A36 = changebasis(P3, P6)
        A67 = changebasis(P6, P7)
        A37 = changebasis(P3, P7)
        A37 ≈ A36 * A67

        @test P2 ⊈ P3

        @test isnondegenerate(P1)
        @test isnondegenerate(P2)
        @test isnondegenerate(P3)
        @test isdegenerate(P4)
        @test changebasis_R(P1, P1) == I
        @test changebasis_R(P2, P2) == I
        @test changebasis_R(P3, P3) == I
        @test changebasis_R(P4, P4) ≠ I
    end

    @testset "sqsubseteq (I)" begin
        p1 = 1
        p2 = 2
        P1 = BSplineSpace{p1}(KnotVector([1, 2, 3, 4, 5]))
        P2 = BSplineSpace{p2}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
        @test P1 ⊑ P1
        @test P1 ⊑ P2
        @test P2 ⊒ P1
        @test P2 ⋢ P1
        @test P1 ⋣ P2

        test_changebasis_I(P1, P2)
        test_changebasis_I(P1, P1)
        test_changebasis_I(P2, P2)

        P3 = BSplineSpace{p1-1}(knotvector(P1)[2:end-1])
        P4 = BSplineSpace{p2-1}(knotvector(P2)[2:end-1])
        @test P3 ⊑ P4
        @test P3 ⊈ P4

        test_changebasis_I(P3, P4)

        # https://github.com/hyrodium/BasicBSpline.jl/issues/325
        P5 = BSplineSpace{2}(knotvector" 21 3")
        P6 = BSplineSpace{3}(knotvector"212 4")
        test_changebasis_I(P5, P6)
    end

    @testset "changebasis_sim" begin
        for p in 1:3, L in 1:8
            k1 = UniformKnotVector(1:L+2p+1)
            k2 = UniformKnotVector(p+1:p+L+1) + p*KnotVector([p+1]) + KnotVector(p+L+1:2p+L)
            P1 = BSplineSpace{p}(k1)
            P2 = BSplineSpace{p}(k2)
            @test P1 isa UniformBSplineSpace{p,T} where {p,T}
            @test domain(P1) == domain(P2)
            @test dim(P1) == dim(P2)
            @test P1 ≃ P2
            A = @inferred BasicBSpline._changebasis_sim(P1,P2)
            D = domain(P1)
            n = dim(P1)
            for _ in 1:100
                t = rand(D)
                @test norm(bsplinebasis.(P1,1:n,t) - A*bsplinebasis.(P2,1:n,t), Inf) < ε
            end
        end
    end

    @testset "non-float" begin
        k = UniformKnotVector(1:15)
        P = BSplineSpace{3}(k)
        @test changebasis(P,P) isa SparseMatrixCSC{Float64}
        @test BasicBSpline._changebasis_R(P,P) isa SparseMatrixCSC{Float64}
        @test BasicBSpline._changebasis_I(P,P) isa SparseMatrixCSC{Float64}
        @test BasicBSpline._changebasis_sim(P,P) isa SparseMatrixCSC{Float64}

        k = UniformKnotVector(1:15//1)
        P = BSplineSpace{3}(k)
        @test changebasis(P,P) isa SparseMatrixCSC{Rational{Int}}
        @test BasicBSpline._changebasis_R(P,P) isa SparseMatrixCSC{Rational{Int}}
        @test BasicBSpline._changebasis_I(P,P) isa SparseMatrixCSC{Rational{Int}}
        @test BasicBSpline._changebasis_sim(P,P) isa SparseMatrixCSC{Rational{Int}}

        k = UniformKnotVector(1:BigInt(15))
        P = BSplineSpace{3}(k)
        @test changebasis(P,P) isa SparseMatrixCSC{BigFloat}
        @test BasicBSpline._changebasis_R(P,P) isa SparseMatrixCSC{BigFloat}
        @test BasicBSpline._changebasis_I(P,P) isa SparseMatrixCSC{BigFloat}
        @test BasicBSpline._changebasis_sim(P,P) isa SparseMatrixCSC{BigFloat}
    end

    @testset "uniform" begin
        for r in 1:5
            k = UniformKnotVector(0:r:50)
            for p in 0:4
                P = BSplineSpace{p}(k)
                Q = BSplineSpace{p,Int,KnotVector{Int}}(P)

                k′ = UniformKnotVector(0:50)
                P′ = BSplineSpace{p}(k′)
                Q′ = BSplineSpace{p,Int,KnotVector{Int}}(P′)
                A1 = @inferred changebasis(P, P′)
                A2 = @inferred changebasis(Q, Q′)
                @test P ⊆ P′
                @test P == Q
                @test P′ == Q′
                @test A1 ≈ A2

                left = leftendpoint(domain(P))-p
                right = rightendpoint(domain(P))+p
                k′ = UniformKnotVector(left:right)
                P′ = BSplineSpace{p}(k′)
                Q′ = BSplineSpace{p,Int,KnotVector{Int}}(P′)
                A3 = @inferred changebasis(P, P′)
                A4 = @inferred changebasis(Q, Q′)
                @test P ⊑ P′
                @test P == Q
                @test P′ == Q′
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
