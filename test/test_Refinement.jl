@testset "Refinement" begin
    p1,p2,p3 = 3,2,4

    k1 = KnotVector(1:8)
    k2 = KnotVector([1,3]) + p2*KnotVector([1,3])
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

    @testset "_I and _R" begin
        @testset "1dim" begin
            a1 = rand(SVector{3, Float64}, n1)
            a2 = rand(SVector{3, Float64}, n2)
            a3 = rand(SVector{3, Float64}, n3)
            C1 = BSplineManifold(a1, P1)
            C2 = BSplineManifold(a2, P2)
            C3 = BSplineManifold(a3, P3)
            @test refinement_R(C1, P1′) == refinement(C1, P1′)
            @test refinement_R(C2, P2′) == refinement(C2, P2′)
            @test_throws DomainError refinement_R(C3, P3′)
            @test_throws DomainError refinement_I(C1, P1′)
            @test refinement_I(C2, P2′) == refinement(C2, P2′)
            @test refinement_I(C3, P3′) == refinement(C3, P3′)
            @test refinement_R(C1) == refinement_I(C1) == refinement(C1)
            @test refinement_R(C2) == refinement_I(C2) == refinement(C2)
            @test refinement_R(C3) == refinement_I(C3) == refinement(C3)
        end

        @testset "2dim" begin
            a12 = rand(SVector{3, Float64}, n1, n2)
            a23 = rand(SVector{3, Float64}, n2, n3)
            S12 = BSplineManifold(a12, P1, P2)
            S23 = BSplineManifold(a23, P2, P3)
            @test refinement_R(S12, P1′, P2′) == refinement(S12, P1′, P2′)
            @test_throws DomainError refinement_R(S23, P2′, P3′) == refinement(S23, P2′, P3′)
            @test_throws DomainError refinement_I(S12, P1′, P2′) == refinement(S12, P1′, P2′)
            @test refinement_I(S23, P2′, P3′) == refinement(S23, P2′, P3′)
        end

        @testset "4dim" begin
            a1122 = rand(SVector{3, Float64}, n1, n1, n2, n2)
            a2233 = rand(SVector{3, Float64}, n2, n2, n3, n3)
            S1122 = BSplineManifold(a1122, P1, P1, P2, P2)
            S2233 = BSplineManifold(a2233, P2, P2, P3, P3)
            @test refinement_R(S1122, P1′, P1′, P2′, P2′) == refinement(S1122, P1′, P1′, P2′, P2′)
            @test_throws DomainError refinement_R(S2233, P2′, P2′, P3′, P3′) == refinement(S2233, P2′, P2′, P3′, P3′)
            @test_throws DomainError refinement_I(S1122, P1′, P1′, P2′, P2′) == refinement(S12, P1′, P1′, P2′, P2′)
            @test refinement_I(S2233, P2′, P2′, P3′, P3′) == refinement(S2233, P2′, P2′, P3′, P3′)
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
        k₊ = (KnotVector([4.4,4.7]),KnotVector(Float64[]),KnotVector([42]))

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
