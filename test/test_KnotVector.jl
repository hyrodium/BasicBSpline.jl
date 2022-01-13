@testset "KnotVector" begin
    k1 = KnotVector([1,2,3])
    k2 = KnotVector([1,2,2,3])
    k3 = KnotVector([2,4,5])
    @testset "constructor" begin
        @test k1 isa KnotVector{Float64}
        @test k1 == KnotVector(1:3)::KnotVector{Float64}
        @test k1 == KnotVector([1,3,2])::KnotVector{Float64}
        @test k1 == KnotVector(1,2,3)::KnotVector{Float64}
        @test k1 == KnotVector(1,3,2)::KnotVector{Float64}
        @test k1 == KnotVector(1.,3,2)::KnotVector{Float64}
        @test k1 == KnotVector{Int}(1,3,2)::KnotVector{Int}
        @test k1 != k2

        @test KnotVector{Int}([1,2]) isa KnotVector{Int}
        @test KnotVector{Int}(1,2) isa KnotVector{Int}
        @test KnotVector{Int}(1,2.) isa KnotVector{Int}
        @test KnotVector{Int}(k1) isa KnotVector{Int}
        @test KnotVector(k1) isa KnotVector{Float64}
    end

    @testset "zeros" begin
        @test KnotVector() == zero(KnotVector)
        @test KnotVector() == 0*k1 == k1*0 == zero(k1)
        @test KnotVector() == KnotVector(Float64[])
    end

    @testset "length" begin
        @test length(KnotVector()) == 0
        @test length(KnotVector([1,2,2,3])) == 4
    end

    @testset "addition, multiply" begin
        @test KnotVector([-1,2,3]) + 2 * KnotVector([2,5]) == KnotVector([-1,2,2,2,3,5,5])
        @test KnotVector([-1,2,3]) + KnotVector([2,5]) * 2 == KnotVector([-1,2,2,2,3,5,5])
        @test k1 + k3 == KnotVector([1,2,2,3,4,5])
        @test 2 * KnotVector([2,3]) == KnotVector([2,2,3,3])

        # type promotion
        @test KnotVector{Int}(1,2) + KnotVector{Rational{Int}}(3) == KnotVector(1,2,3)
        @test KnotVector{Int}(1,2) + KnotVector{Rational{Int}}(3) isa KnotVector{Rational{Int}}
        @test KnotVector{Int}(1,2)*0 == KnotVector()
        @test KnotVector{Int}(1,2)*0 isa KnotVector{Int}
        @test KnotVector{Int}() isa KnotVector{Int}
    end

    @testset "unique" begin
        @test unique(k1) == k1
        @test unique(k2) == k1
    end

    @testset "inclusive relation" begin
        @test KnotVector([1,2,3])     ⊆ KnotVector([1,2,3])
        @test KnotVector([1,2,3])     ⊇ KnotVector([1,2,3])
        @test KnotVector([1,2,2,3])   ⊆ KnotVector([1,2,2,3,5])
        @test KnotVector([1,2,2,3])   ⊈ KnotVector([1,2,3,5])
        @test KnotVector([1,2,2,3,5]) ⊇ KnotVector([1,2,2,3])
        @test KnotVector([1,2,3,5])   ⊉ KnotVector([1,2,2,3])
    end

    @testset "string" begin
        k = KnotVector([1,2,2,3])
        @test string(k) == "KnotVector([1.0, 2.0, 2.0, 3.0])"
        @test string(KnotVector()) == "KnotVector([])"
    end

    @testset "other operators" begin
        k = KnotVector([1,2,2,3])
        @test 𝔫(k, 0.3) == 0
        @test 𝔫(k, 1.0) == 1
        @test 𝔫(k, 2.0) == 2
        @test 1 ∈ k
        @test 1.5 ∉ k
    end
end
