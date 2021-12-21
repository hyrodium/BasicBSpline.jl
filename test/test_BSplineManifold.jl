@testset "BSplineManifold" begin
    function arrayofvector2array(a::AbstractArray{<:AbstractVector{R}, d})::Array{R,d+1} where {d, R<:Real}
        d̂ = length(a[1])
        s = size(a)
        N = prod(s)
        a_2dim = [a[i][j] for i in 1:N, j in 1:d̂]
        a′ = reshape(a_2dim, s..., d̂)
        return a′
    end

    function array2arrayofvector(a::AbstractArray{R,d})::Array{Vector{R}, d-1} where {d, R<:Real}
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

            P1 = BSplineSpace{1}(Knots([0, 0, 1, 1]))
            n1 = dim(P1) # 2
            ax = [i for i in 1:n1]
            ay = [rand() for i in 1:n1]
            a = hcat(ax,ay)  # n1 × 2 array of d̂-dim vector.
            M = BSplineManifold((P1,), a)
            @test dim(M) == 1

            P1′ = BSplineSpace{2}(Knots([-2, 0, 0, 1, 1, 2]))
            p₊ = (1,)
            k₊ = (Knots(0.3),)

            @test P1 ⊑ P1′

            M′ = refinement(M, (P1′,))
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [rand() for _ in 1:10]
            for t in ts
                @test M(t) ≈ M′(t)
                @test M(t) ≈ M′′(t)
            end
        end

        @testset "FastBSplineManifold-1dim" begin
            Random.seed!(42)

            P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
            n1 = dim(P1) # 2
            a = [Point(i, rand()) for i in 1:n1]  # n1 × n2 array of d̂-dim vector.
            M = FastBSplineManifold([P1], a)
            @test dim(M) == 1

            P1′ = FastBSplineSpace(2, Knots([-2, 0, 0, 1, 1, 2]))
            p₊ = [1]
            k₊ = [Knots()]

            @test P1 ⊑ P1′

            M′ = refinement(M, [P1′])
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [[rand()] for _ in 1:10]
            for t in ts
                @test M(t...) ≈ M′(t...)
                @test M(t...) ≈ M′′(t...)
            end
        end

        @testset "BSplineCurve" begin
            Random.seed!(42)

            P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
            n1 = dim(P1) # 2
            a = [Point(i, rand()) for i in 1:n1]  # n1 × n2 array of d̂-dim vector.
            M = BSplineCurve([P1], a)
            @test dim(M) == 1

            P1′ = FastBSplineSpace(2, Knots([-2, 0, 0, 1, 1, 2]))
            p₊ = [1]
            k₊ = [Knots()]

            @test P1 ⊑ P1′

            M′ = refinement(M, [P1′])
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [[rand()] for _ in 1:10]
            for t in ts
                @test M(t) ≈ M′(t)
                @test M(t) ≈ M′′(t)
            end
        end
    end

    @testset "2dim" begin
        @testset "BSplineManifold-2dim" begin
            Random.seed!(42)

            P1 = BSplineSpace{1}(Knots([0, 0, 1, 1]))
            P2 = BSplineSpace{1}(Knots([1, 1, 2, 3, 3]))
            n1 = dim(P1) # 2
            n2 = dim(P2) # 3
            rand_a = [rand(2) for i in 1:n1, j in 1:n2]
            a = [[2*i-6.5,2*j-6.5] for i in 1:n1, j in 1:n2] + rand_a
            M = BSplineManifold((P1, P2), arrayofvector2array(a))
            @test dim(M) == 2

            P1′ = BSplineSpace{2}(Knots([-2, 0, 0, 1, 1, 2]))
            P2′ = BSplineSpace{1}(Knots([-3, 1, 2, 1.45, 3, 4]))
            p₊ = (1, 0)
            k₊ = (Knots(), Knots(1.45))

            @test P1 ⊑ P1′
            @test P2 ⊑ P2′

            M′ = refinement(M, (P1′, P2′))
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
            for t in ts
                @test M(t...) ≈ M′(t...)
                @test M(t...) ≈ M′′(t...)
            end
        end

        @testset "FastBSplineManifold-2dim" begin
            Random.seed!(42)

            P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
            P2 = FastBSplineSpace(1, Knots([1, 1, 2, 3, 3]))
            n1 = dim(P1) # 2
            n2 = dim(P2) # 3
            a = [Point(i, j) for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
            M = FastBSplineManifold([P1, P2], a)
            @test dim(M) == 2

            P1′ = FastBSplineSpace(2, Knots([0, 0, 0, 1, 1, 1]))
            P2′ = FastBSplineSpace(1, Knots([1, 1, 2, 1.45, 3, 3]))
            p₊ = [1, 0]
            k₊ = [Knots(), Knots(1.45)]

            @test P1 ⊆ P1′
            @test P2 ⊆ P2′

            M′ = refinement(M, [P1′, P2′])
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
            for t in ts
                @test M(t...) ≈ M′(t...)
                @test M(t...) ≈ M′′(t...)
            end
        end

        @testset "BSplineSurface" begin
            Random.seed!(42)

            P1 = FastBSplineSpace(1, Knots([0, 0, 1, 1]))
            P2 = FastBSplineSpace(1, Knots([1, 1, 2, 3, 3]))
            n1 = dim(P1) # 2
            n2 = dim(P2) # 3
            a = [Point(i, j) for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂-dim vector.
            M = BSplineSurface([P1, P2], a)
            @test dim(M) == 2

            P1′ = FastBSplineSpace(2, Knots([0, 0, 0, 1, 1, 1]))
            P2′ = FastBSplineSpace(1, Knots([1, 1, 2, 1.45, 3, 3]))
            p₊ = [1, 0]
            k₊ = [Knots(), Knots(1.45)]

            @test P1 ⊆ P1′
            @test P2 ⊆ P2′

            M′ = refinement(M, [P1′, P2′])
            M′′ = refinement(M, p₊=p₊, k₊=k₊)
            ts = [[rand(), 1 + 2 * rand()] for _ in 1:10]
            for t in ts
                @test M(t) ≈ M′(t)
                @test M(t) ≈ M′′(t)
            end
        end
    end
end
