@testset "Refinement" begin
    p1,p2,p3 = 3,2,4

    k1 = KnotVector(1:8)
    k2 = KnotVector([1,4]) + p2*KnotVector([1,4])
    k3 = KnotVector([0,2]) + p3*KnotVector([1,2])

    P1 = BSplineSpace{p1}(k1)
    P2 = BSplineSpace{p2}(k2)
    P3 = BSplineSpace{p3}(k3)

    P1′ = BSplineSpace{3}(k1+KnotVector([1.2,5.5,7.2,8.5]))
    P2′ = expandspace(P2, Val(1), KnotVector([2,2,2,2,2]))
    P3′ = BSplineSpace{4}(UniformKnotVector(-3.0:6.0))

    n1,n2,n3 = dim(P1),dim(P2),dim(P3)
    D1,D2,D3 = domain(P1),domain(P2),domain(P3)

    @test P1 ⊆ P1′
    @test P2 ⊆ P2′
    @test !(P3 ⊆ P3′)

    @test !(P1 ⊑ P1′)
    @test P2 ⊑ P2′
    @test P3 ⊑ P3′

    @test P3 ≃ P3′

    @testset "_I and _R" begin
        @testset "1dim" begin
            a1 = rand(SVector{3, Float64}, n1)
            a2 = rand(SVector{3, Float64}, n2)
            a3 = rand(SVector{3, Float64}, n3)
            w1 = rand(n1) .+ 1
            w2 = rand(n2) .+ 1
            w3 = rand(n3) .+ 1
            M1 = BSplineManifold(a1, P1)
            M2 = BSplineManifold(a2, P2)
            M3 = BSplineManifold(a3, P3)
            R1 = RationalBSplineManifold(a1, w1, P1)
            R2 = RationalBSplineManifold(a2, w2, P2)
            R3 = RationalBSplineManifold(a3, w3, P3)

            @test refinement_R(M1, P1′) == refinement(M1, P1′)
            @test refinement_R(R1, P1′) == refinement(R1, P1′)
            @test refinement_R(M2, P2′) == refinement(M2, P2′)
            @test refinement_R(R2, P2′) == refinement(R2, P2′)
            @test_throws DomainError refinement_R(M3, P3′)
            @test_throws DomainError refinement_R(R3, P3′)
            @test_throws DomainError refinement_I(M1, P1′)
            @test_throws DomainError refinement_I(R1, P1′)
            @test refinement_I(M2, P2′) == refinement(M2, P2′)
            @test refinement_I(R2, P2′) == refinement(R2, P2′)
            @test refinement_I(M3, P3′) == refinement(M3, P3′)
            @test refinement_I(R3, P3′) == refinement(R3, P3′)

            @test refinement_R(M1) == refinement_I(M1) == refinement(M1)
            @test refinement_R(R1) == refinement_I(R1) == refinement(R1)
            @test refinement_R(M2) == refinement_I(M2) == refinement(M2)
            @test refinement_R(R2) == refinement_I(R2) == refinement(R2)
            @test refinement_R(M3) == refinement_I(M3) == refinement(M3)
            @test refinement_R(R3) == refinement_I(R3) == refinement(R3)

            # P1 ⊆ P1′ but not P1 ⊑ P1′
            M1a_R = refinement_R(refinement_R(M1, (Val(0),)), (KnotVector([1.2,5.5,7.2,8.5]),))
            M1b_R = refinement_R(M1, (Val(0),), (KnotVector([1.2,5.5,7.2,8.5]),))
            M1c_R = refinement_R(M1, P1′)
            @test controlpoints(M1a_R) ≈ controlpoints(M1b_R) == controlpoints(M1c_R)
            @test_throws DomainError refinement_I(refinement_I(M1, (Val(0),)), (KnotVector([1.2,5.5,7.2,8.5]),))
            @test_throws DomainError refinement_I(M1, (Val(0),), (KnotVector([1.2,5.5,7.2,8.5]),))
            @test_throws DomainError refinement_I(M1, P1′)
            @test_throws DomainError refinement(refinement(M1, (Val(1),)), (KnotVector([1.2,5.5,7.2,8.5]),))
            @test_throws DomainError refinement(M1, (Val(1),), (KnotVector([1.2,5.5,7.2,8.5]),))
            M1c = refinement(M1, P1′)

            # P2 ⊆ P2′ and P2 ⊑ P2′
            M2a_R = refinement_R(refinement_R(M2, (Val(1),)), (KnotVector([2,2,2,2,2]),))
            M2b_R = refinement_R(M2, (Val(1),), (KnotVector([2,2,2,2,2]),))
            M2c_R = refinement_R(M2, P2′)
            @test controlpoints(M2a_R) ≈ controlpoints(M2b_R) ≠ controlpoints(M2c_R)
            M2a_I = refinement_I(refinement_I(M2, (Val(1),)), (KnotVector([2,2,2,2,2]),))
            M2b_I = refinement_I(M2, (Val(1),), (KnotVector([2,2,2,2,2]),))
            M2c_I = refinement_I(M2, P2′)
            @test controlpoints(M2a_I) ≈ controlpoints(M2b_I) == controlpoints(M2c_I)
            M2a = refinement(refinement(M2, (Val(1),)), (KnotVector([2,2,2,2,2]),))
            M2b = refinement(M2, (Val(1),), (KnotVector([2,2,2,2,2]),))
            M2c = refinement(M2, P2′)
            @test controlpoints(M2a) ≈ controlpoints(M2b) == controlpoints(M2c)

            # P3 ⊑ P3′ but not P3 ⊆ P3′
            M3a_R = refinement_R(refinement_R(M3, (Val(0),)), (EmptyKnotVector(),))
            M3b_R = refinement_R(M3, (Val(0),), (EmptyKnotVector(),))
            @test_throws DomainError refinement_R(M3, P3′)
            @test controlpoints(M3a_R) ≈ controlpoints(M3b_R)
            M3a_I = refinement_I(refinement_I(M3, (Val(0),)), (EmptyKnotVector(),))
            M3b_I = refinement_I(M3, (Val(0),), (EmptyKnotVector(),))
            M3c_I = refinement_I(M3, P3′)
            @test controlpoints(M3a_I) ≈ controlpoints(M3b_I) ≉ controlpoints(M3c_I)
            M3a = refinement(refinement(M3, (Val(0),)), (EmptyKnotVector(),))
            M3b = refinement(M3, (Val(0),), (EmptyKnotVector(),))
            M3c = refinement(M3, P3′)
            @test controlpoints(M3a) ≈ controlpoints(M3b) ≉ controlpoints(M3c)
        end

        @testset "2dim" begin
            a12 = rand(SVector{3, Float64}, n1, n2)
            a23 = rand(SVector{3, Float64}, n2, n3)
            w12 = rand(n1, n2) .+ 1
            w23 = rand(n2, n3) .+ 1
            M12 = BSplineManifold(a12, P1, P2)
            M23 = BSplineManifold(a23, P2, P3)
            R12 = RationalBSplineManifold(a12, w12, P1, P2)
            R23 = RationalBSplineManifold(a23, w23, P2, P3)

            @test refinement_R(M12, P1′, P2′) == refinement(M12, P1′, P2′)
            @test refinement_R(R12, P1′, P2′) == refinement(R12, P1′, P2′)
            @test_throws DomainError refinement_R(M23, P2′, P3′) == refinement(M23, P2′, P3′)
            @test_throws DomainError refinement_R(R23, P2′, P3′) == refinement(R23, P2′, P3′)
            @test_throws DomainError refinement_I(M12, P1′, P2′) == refinement(M12, P1′, P2′)
            @test_throws DomainError refinement_I(R12, P1′, P2′) == refinement(R12, P1′, P2′)
            @test refinement_I(M23, P2′, P3′) == refinement(M23, P2′, P3′)
            @test refinement_I(R23, P2′, P3′) == refinement(R23, P2′, P3′)

            @test refinement_R(M12) == refinement_I(M12) == refinement(M12)
            @test refinement_R(R12) == refinement_I(R12) == refinement(R12)
            @test refinement_R(M23) == refinement_I(M23) == refinement(M23)
            @test refinement_R(R23) == refinement_I(R23) == refinement(R23)
        end

        @testset "4dim" begin
            a1122 = rand(SVector{3, Float64}, n1, n1, n2, n2)
            a2233 = rand(SVector{3, Float64}, n2, n2, n3, n3)
            w1122 = rand(n1, n1, n2, n2) .+ 1
            w2233 = rand(n2, n2, n3, n3) .+ 1
            M1122 = BSplineManifold(a1122, P1, P1, P2, P2)
            M2233 = BSplineManifold(a2233, P2, P2, P3, P3)
            R1122 = RationalBSplineManifold(a1122, w1122, P1, P1, P2, P2)
            R2233 = RationalBSplineManifold(a2233, w2233, P2, P2, P3, P3)

            @test refinement_R(M1122, P1′, P1′, P2′, P2′) == refinement(M1122, P1′, P1′, P2′, P2′)
            @test refinement_R(R1122, P1′, P1′, P2′, P2′) == refinement(R1122, P1′, P1′, P2′, P2′)
            @test_throws DomainError refinement_R(M2233, P2′, P2′, P3′, P3′) == refinement(M2233, P2′, P2′, P3′, P3′)
            @test_throws DomainError refinement_R(R2233, P2′, P2′, P3′, P3′) == refinement(R2233, P2′, P2′, P3′, P3′)
            @test_throws DomainError refinement_I(M1122, P1′, P1′, P2′, P2′) == refinement(M12, P1′, P1′, P2′, P2′)
            @test_throws DomainError refinement_I(R1122, P1′, P1′, P2′, P2′) == refinement(R12, P1′, P1′, P2′, P2′)
            @test refinement_I(M2233, P2′, P2′, P3′, P3′) == refinement(M2233, P2′, P2′, P3′, P3′)
            @test refinement_I(R2233, P2′, P2′, P3′, P3′) == refinement(R2233, P2′, P2′, P3′, P3′)

            @test refinement_R(M1122) == refinement_I(M1122) == refinement(M1122)
            @test refinement_R(R1122) == refinement_I(R1122) == refinement(R1122)
            @test refinement_R(M2233) == refinement_I(M2233) == refinement(M2233)
            @test refinement_R(R2233) == refinement_I(R2233) == refinement(R2233)
        end
    end

    @testset "1dim" begin
        a = [Point(i, rand()) for i in 1:n1]
        p₊ = (Val(1),)
        k₊ = (KnotVector([4.5]),)
        # k₊ = (KnotVector([4.5,4.95]),)

        @testset "BSplineManifold" begin
            M = BSplineManifold(a, (P1,))
            M0 = @inferred refinement(M)
            M1 = @inferred refinement(M, k₊)
            M2 = @inferred refinement(M, p₊)
            M3 = @inferred refinement(M, p₊, k₊)
            M4 = @inferred refinement(M, (P1′,))

            for _ in 1:100
                t1 = rand(D1)
                @test M(t1) ≈ M1(t1)
                @test M(t1) ≈ M2(t1)
                @test M(t1) ≈ M3(t1)
                @test M(t1) ≈ M4(t1)
            end
        end

        @testset "RationalBSplineManifold" begin
            w = rand(n1) .+ 1
            R = RationalBSplineManifold(a, w, (P1,))
            R0 = @inferred refinement(R)
            R1 = @inferred refinement(R, k₊)
            R2 = @inferred refinement(R, p₊)
            R3 = @inferred refinement(R, p₊, k₊)
            R4 = @inferred refinement(R, (P1′,))

            for _ in 1:100
                t1 = rand(D1)
                @test R(t1) ≈ R1(t1)
                @test R(t1) ≈ R2(t1)
                @test R(t1) ≈ R3(t1)
                @test R(t1) ≈ R4(t1)
            end
        end
    end

    @testset "2dim" begin
        a = [Point(i, rand()) for i in 1:n1, j in 1:n2]
        p₊ = (Val(1), Val(2))
        k₊ = (KnotVector([4.5,4.7]),KnotVector(Float64[]))

        @testset "BSplineManifold" begin
            M = BSplineManifold(a, (P1,P2))
            M0 = @inferred refinement(M)
            M1 = @inferred refinement(M, k₊)
            M2 = @inferred refinement(M, p₊)
            M3 = @inferred refinement(M, p₊, k₊)
            M4 = @inferred refinement(M, (P1′,P2′))

            for _ in 1:100
                t1 = rand(D1)
                t2 = rand(D2)
                @test M(t1,t2) ≈ M1(t1,t2)
                @test M(t1,t2) ≈ M2(t1,t2)
                @test M(t1,t2) ≈ M3(t1,t2)
                @test M(t1,t2) ≈ M4(t1,t2)
            end
        end

        @testset "RationalBSplineManifold" begin
            w = rand(n1,n2) .+ 1

            R = RationalBSplineManifold(a, w, (P1,P2))
            R0 = @inferred refinement(R)
            R1 = @inferred refinement(R, k₊)
            R2 = @inferred refinement(R, p₊)
            R3 = @inferred refinement(R, p₊, k₊)
            R4 = @inferred refinement(R, (P1′,P2′))

            for _ in 1:100
                t1 = rand(D1)
                t2 = rand(D2)
                @test R(t1,t2) ≈ R1(t1,t2)
                @test R(t1,t2) ≈ R2(t1,t2)
                @test R(t1,t2) ≈ R3(t1,t2)
                @test R(t1,t2) ≈ R4(t1,t2)
            end
        end
    end

    @testset "3dim" begin
        a = [Point(i, rand()) for i in 1:n1, j in 1:n2, k in 1:n3]
        p₊ = (Val(1), Val(2), Val(0))
        k₊ = (KnotVector([4.4,4.7]),KnotVector(Float64[]),KnotVector([1.8]))

        @testset "BSplineManifold" begin
            M = BSplineManifold(a, (P1,P2,P3))
            # On Julia v1.6, the following script seems not type-stable.
            if VERSION ≥ v"1.8"
                M0 = @inferred refinement(M)
                M1 = @inferred refinement(M, k₊)
                M2 = @inferred refinement(M, p₊)
                M3 = @inferred refinement(M, p₊, k₊)
                M4 = @inferred refinement(M, (P1′,P2′,P3′))
            else
                M0 = refinement(M)
                M1 = refinement(M, k₊)
                M2 = refinement(M, p₊)
                M3 = refinement(M, p₊, k₊)
                M4 = refinement(M, (P1′,P2′,P3′))
            end

            for _ in 1:100
                t1 = rand(D1)
                t2 = rand(D2)
                t3 = rand(D3)
                @test M(t1,t2,t3) ≈ M1(t1,t2,t3)
                @test M(t1,t2,t3) ≈ M2(t1,t2,t3)
                @test M(t1,t2,t3) ≈ M3(t1,t2,t3)
                @test M(t1,t2,t3) ≈ M4(t1,t2,t3)
            end
        end

        @testset "RationalBSplineManifold" begin
            w = rand(n1,n2,n3) .+ 1

            R = RationalBSplineManifold(a, w, (P1,P2,P3))
            # On Julia v1.6, the following script seems not type-stable.
            if VERSION ≥ v"1.8"
                R0 = @inferred refinement(R)
                R1 = @inferred refinement(R, k₊)
                R2 = @inferred refinement(R, p₊)
                R3 = @inferred refinement(R, p₊, k₊)
                R4 = @inferred refinement(R, (P1′,P2′,P3′))
            else
                R0 = refinement(R)
                R1 = refinement(R, k₊)
                R2 = refinement(R, p₊)
                R3 = refinement(R, p₊, k₊)
                R4 = refinement(R, (P1′,P2′,P3′))
            end

            for _ in 1:100
                t1 = rand(D1)
                t2 = rand(D2)
                t3 = rand(D3)
                @test R(t1,t2,t3) ≈ R1(t1,t2,t3)
                @test R(t1,t2,t3) ≈ R2(t1,t2,t3)
                @test R(t1,t2,t3) ≈ R3(t1,t2,t3)
                @test R(t1,t2,t3) ≈ R4(t1,t2,t3)
            end
        end
    end
end
