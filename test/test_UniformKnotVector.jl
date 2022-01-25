@testset "UniformKnotVector" begin
    k1 = UniformKnotVector(1:1:3)
    k2 = UniformKnotVector(1:3)
    k3 = UniformKnotVector(1:4)
    k4 = UniformKnotVector(2:4)

    @testset "equality" begin
        @test k1 == k1
        @test k1 == k2
        @test k2 != k3
        @test k4 != k3
        @test k4 != k2
    end

    @testset "constructor" begin
        @test k2 == KnotVector(1:3)::KnotVector{Float64}
        @test k2 isa UniformKnotVector{Int,UnitRange{Int}}
        @test UniformKnotVector(k2) isa UniformKnotVector{Int,UnitRange{Int}}
        @test_throws MethodError UniformKnotVector{Float64}(k3)
        T = UniformKnotVector{Float64,StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
        @test T(k3) isa T
    end

    @testset "length" begin
        @test length(k2) == 3
        @test length(k3) == 4
        @test length(k4) == 3
    end

    @testset "addition, multiply" begin
        @test k2 + k3 == KnotVector(k2) + KnotVector(k3)
        @test k2 + k4 == KnotVector(k2) + KnotVector(k4)
        @test 2 * k2 == k2 * 2 == k2 + k2
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

    # @testset "string" begin
    #     k = UniformKnotVector([1,2,2,3])
    #     @test string(k) == "UniformKnotVector([1.0, 2.0, 2.0, 3.0])"
    #     @test string(UniformKnotVector()) == "UniformKnotVector([])"
    # end

    @testset "other operators" begin
        @test 𝔫(k2, 0.3) == 0
        @test 𝔫(k2, 1) == 1
        @test 𝔫(k2, 1.0) == 1
        @test 𝔫(k2, 2.0) == 1
        @test 1 ∈ k2
        @test 1.5 ∉ k2
    end
end
