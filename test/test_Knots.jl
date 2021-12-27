@testset "Knots" begin
    k1 = Knots([1,2,3])
    k2 = Knots([1,2,2,3])
    k3 = Knots([2,4,5])
    @testset "constructor" begin
        @test k1 isa Knots{Float64}
        @test k1 == Knots(1:3)::Knots{Float64}
        @test k1 == Knots([1,3,2])::Knots{Float64}
        @test k1 == Knots(1,2,3)::Knots{Float64}
        @test k1 == Knots(1,3,2)::Knots{Float64}

        @test Knots{Int}([1,2]) isa Knots{Int}
        @test Knots{Int}(1,2) isa Knots{Int}
    end

    @testset "zeros" begin
        @test Knots() == zero(Knots)
        @test Knots() == 0*k1
        @test Knots() == Knots(Float64[])
    end

    @testset "length" begin
        @test length(Knots()) == 0
        @test length(Knots([1,2,2,3])) == 4
    end

    @testset "addition, multiply" begin
        @test Knots([-1,2,3]) + 2 * Knots([2,5]) == Knots([-1,2,2,2,3,5,5])
        @test k1 + k3 == Knots([1,2,2,3,4,5])
        @test 2 * Knots([2,3]) == Knots([2,2,3,3])
    end

    @testset "unique" begin
        @test unique(k1) == k1
        @test unique(k2) == k1
    end

    @testset "inclusive relation" begin
        @test Knots([1,2,3])     ⊆ Knots([1,2,3])
        @test Knots([1,2,3])     ⊇ Knots([1,2,3])
        @test Knots([1,2,2,3])   ⊆ Knots([1,2,2,3,5])
        @test Knots([1,2,2,3])   ⊈ Knots([1,2,3,5])
        @test Knots([1,2,2,3,5]) ⊇ Knots([1,2,2,3])
        @test Knots([1,2,3,5])   ⊉ Knots([1,2,2,3])
    end

    @testset "string" begin
        k = Knots([1,2,2,3])
        @test string(k) == "Knots([1.0, 2.0, 2.0, 3.0])"
        @test string(Knots()) == "Knots([])"
    end

    @testset "other operators" begin
        k = Knots([1,2,2,3])
        @test 𝔫(k, 0.3) == 0
        @test 𝔫(k, 1.0) == 1
        @test 𝔫(k, 2.0) == 2
        @test 1 ∈ k
        @test 1.5 ∉ k
    end
end
