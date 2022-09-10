@testset "BSplineManifold" begin
    function arrayofvector2array(a::AbstractArray{<:AbstractVector{T},d})::Array{T,d+1} where {d, T<:Real}
        d̂ = length(a[1])
        s = size(a)
        N = prod(s)
        a_2dim = [a[i][j] for i in 1:N, j in 1:d̂]
        a′ = reshape(a_2dim, s..., d̂)
        return a′
    end

    function array2arrayofvector(a::AbstractArray{T,d})::Array{Vector{T},d-1} where {d, T<:Real}
        s = size(a)
        d̂ = s[end]
        N = s[1:end-1]
        a_flat = reshape(a,prod(N),d̂)
        a_vec = [a_flat[i,:] for i in 1:prod(N)]
        a′ = reshape(a_vec,N)
        return a′
    end

    @testset "1dim" begin
        @testset "BSplineManifold-1dim" begin
            Random.seed!(42)

            P1 = BSplineSpace{1}(KnotVector([0, 0, 1, 1]))
            n1 = dim(P1) # 2
            a = [Point(i, rand()) for i in 1:n1]  # n1 × n2 array of d̂-dim vector.
            M = BSplineManifold(a, (P1,))
            @test dim(M) == 1

            P1′ = BSplineSpace{2}(KnotVector([-2, 0, 0, 1, 1, 2]))
            p₊ = (Val(1),)
            k₊ = (KnotVector(Float64[]),)

            @test P1 ⊑ P1′

            M′ = refinement(M, P1′)
            M′′ = refinement(M, p₊, k₊)
            ts = [[rand()] for _ in 1:10]
            for t in ts
                @test M(t...) ≈ M′(t...)
                @test M(t...) ≈ M′′(t...)
            end
            @test_throws DomainError M(-5)

            @test map.(M,[0.2,0.5,0.6]) == M.([0.2,0.5,0.6])
            @test M(:) == M
            @test Base.mightalias(controlpoints(M), controlpoints(M))
            @test !Base.mightalias(controlpoints(M), controlpoints(M(:)))
        end
    end

    @testset "2dim" begin
        @testset "BSplineManifold-2dim" begin
            Random.seed!(42)

            P1 = BSplineSpace{1}(KnotVector([0, 0, 1, 1]))
            P2 = BSplineSpace{1}(KnotVector([1, 1, 2, 3, 3]))
            n1 = dim(P1) # 2
            n2 = dim(P2) # 3
            a = [Point(i, j) for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
            M = BSplineManifold(a, (P1, P2))
            @test dim(M) == 2

            P1′ = BSplineSpace{2}(KnotVector([0, 0, 0, 1, 1, 1]))
            P2′ = BSplineSpace{1}(KnotVector([1, 1, 2, 1.45, 3, 3]))
            p₊ = (Val(1), Val(0))
            k₊ = (KnotVector(Float64[]), KnotVector([1.45]))

            @test P1 ⊆ P1′
            @test P2 ⊆ P2′

            M′ = refinement(M, (P1′, P2′))
            M′′ = refinement(M, p₊, k₊)
            ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
            for t in ts
                t1, t2 = t
                @test M(t1,t2) ≈ M′(t1,t2)
                @test M(t1,t2) ≈ M′′(t1,t2)
                @test M(t1,t2) ≈ M(t1,:)(t2)
                @test M(t1,t2) ≈ M(:,t2)(t1)
            end
            @test_throws DomainError M(-5,-8)

            @test M(:,:) == M
            @test Base.mightalias(controlpoints(M), controlpoints(M))
            @test !Base.mightalias(controlpoints(M), controlpoints(M(:,:)))
        end
    end

    @testset "3dim" begin
        @testset "BSplineManifold-3dim" begin
            Random.seed!(42)

            P1 = BSplineSpace{1}(KnotVector([0, 0, 1, 1]))
            P2 = BSplineSpace{1}(KnotVector([1, 1, 2, 3, 3]))
            P3 = BSplineSpace{2}(KnotVector(rand(16)))
            n1 = dim(P1) # 2
            n2 = dim(P2) # 3
            n3 = dim(P3)
            a = [Point(i1, i2, 3) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
            M = BSplineManifold(a, (P1, P2, P3))
            @test dim(M) == 3

            p₊ = (Val(1), Val(0), Val(1))
            k₊ = (KnotVector(Float64[]), KnotVector([1.45]), EmptyKnotVector())

            M′′ = refinement(M, p₊, k₊)
            ts = [[rand(), 1 + 2 * rand(), 0.5] for _ in 1:10]
            for t in ts
                t1, t2, t3 = t
                @test M(t1,t2,t3) ≈ M′′(t1,t2,t3)
                @test M(t1,t2,t3) ≈ M(t1,:,:)(t2,t3)
                @test M(t1,t2,t3) ≈ M(:,t2,:)(t1,t3)
                @test M(t1,t2,t3) ≈ M(:,:,t3)(t1,t2)
                @test M(t1,t2,t3) ≈ M(t1,t2,:)(t3)
                @test M(t1,t2,t3) ≈ M(t1,:,t3)(t2)
                @test M(t1,t2,t3) ≈ M(:,t2,t3)(t1)
            end
            @test_throws DomainError M(-5,-8,-120)

            @test M(:,:,:) == M
            @test Base.mightalias(controlpoints(M), controlpoints(M))
            @test !Base.mightalias(controlpoints(M), controlpoints(M(:,:,:)))
        end
    end
end
