@testset "UniformKnotVector" begin
    k1 = UniformKnotVector(1:1:3)
    k2 = UniformKnotVector(Base.OneTo(3))
    k3 = UniformKnotVector(1:4)
    k4 = UniformKnotVector(2:4)
    k5 = UniformKnotVector(2.0:4.0)
    k6 = UniformKnotVector(2.0:2.0:12.0)

    @testset "equality" begin
        @test k1 == k1
        @test k1 == k2
        @test k2 != k3
        @test k4 != k3
        @test k4 != k2
        @test k2 == KnotVector(1:3)
        @test k2 != 1:3
    end

    @testset "constructor, conversion" begin
        @test k1 isa UniformKnotVector{Int,StepRange{Int,Int}}
        @test k2 isa UniformKnotVector{Int,Base.OneTo{Int}}
        @test k3 isa UniformKnotVector{Int,UnitRange{Int}}
        @test UniformKnotVector(k2) isa UniformKnotVector{Int,Base.OneTo{Int}}
        @test_throws MethodError UniformKnotVector{Float64}(k3)
        T1 = UniformKnotVector{Float64,typeof(1.0:1.0)}
        T2 = UniformKnotVector{Int,typeof(1:1:1)}
        T3 = UniformKnotVector{Float64,typeof(1:1:1)}
        @test T1 != T2
        @test T1(k3) isa T1
        @test T2(k3) isa T2
        @test_throws MethodError T3(k3)
        @test convert(T1,k3) isa T1
        @test convert(T2,k3) isa T2
        @test_throws MethodError convert(T3,k3)
        @test KnotVector(k1) isa KnotVector{Int}
        @test KnotVector(k2) isa KnotVector{Int}
        @test KnotVector(k5) isa KnotVector{Float64}
    end

    @testset "eltype" begin
        @test eltype(UniformKnotVector(1:3)) == Int
        @test eltype(UniformKnotVector(1:1:3)) == Int
        @test eltype(UniformKnotVector(1:1:BigInt(3))) == BigInt
        @test eltype(UniformKnotVector(1.0:1.0:3.0)) == Float64
        @test eltype(UniformKnotVector(1//1:1//5:3//1)) == Rational{Int}
    end

    @testset "length" begin
        @test length(k2) == 3
        @test length(k3) == 4
        @test length(k4) == 3
    end

    @testset "iterator, getindex" begin
        for t in k2
            @test t in k3
        end
        @test k1[1] == 1
        @test k2[end] == 3
        @test k4[2] == 3
        @test k3[2:4] == k4
        @test collect(k2) isa Vector{Int}
        @test [k2...] isa Vector{Int}
        @test collect(k5) isa Vector{Float64}
        @test [k5...] isa Vector{Float64}
        @test collect(k2) == [k2...]
        @test collect(k2) != k2
    end

    @testset "addition, multiply" begin
        @test k2 + k3 == KnotVector(k2) + KnotVector(k3)
        @test k2 + k4 == KnotVector(k2) + KnotVector(k4)
        @test 2 * k2 == k2 * 2 == k2 + k2
    end

    @testset "zeros" begin
        @test_throws MethodError zero(UniformKnotVector)
        @test_throws MethodError zeros(UniformKnotVector,3)
        @test zero(k1) isa KnotVector{Int}
        @test zero(k5) isa KnotVector{Float64}
        @test k1 * 0 == zero(k1)
        @test k1 + zero(k1) == k1
    end

    @testset "unique" begin
        @test unique(k2) == k2
        @test unique(k3) == k3
        @test unique(k4) == k4
    end

    @testset "inclusive relation" begin
        @test k1 == k2
        @test k1 ⊆ k2
        @test k2 ⊆ k1
        @test k2 ⊆ k3
        @test k4 ⊆ k3
        @test k2 ⊈ k4
        @test KnotVector(k2) ⊆ k3
        @test KnotVector(k4) ⊆ k3
        @test KnotVector(k2) ⊈ k4
        @test k2 ⊆ KnotVector(k3)
        @test k4 ⊆ KnotVector(k3)
        @test k2 ⊈ KnotVector(k4)
    end

    @testset "string" begin
        @test string(k1) == "UniformKnotVector(1:1:3)"
        @test string(k2) == "UniformKnotVector(Base.OneTo(3))"
        @test string(k3) == "UniformKnotVector(1:4)"
    end

    @testset "other operators" begin
        @test 𝔫(k2, 0.3) == 0
        @test 𝔫(k2, 1) == 1
        @test 𝔫(k2, 1.0) == 1
        @test 𝔫(k2, 2.0) == 1
        @test 1 ∈ k2
        @test 1.5 ∉ k2
    end
end
