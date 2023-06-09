@testset "KnotVector" begin
    k1 = KnotVector([1,2,3])
    k2 = KnotVector([1,2,2,3])
    k3 = KnotVector([2,4,5])
    @testset "constructor" begin
        @test k1 isa KnotVector{Int64}
        @test k1 == KnotVector(1:3)::KnotVector{Int}
        @test k1 == KnotVector([1,3,2])::KnotVector{Int64}
        @test k1 == KnotVector([1,2,3])::KnotVector{Int64}
        @test k1 == KnotVector([1,3,2])::KnotVector{Int64}
        @test k1 == KnotVector([1.,3,2])::KnotVector{Float64}
        @test k1 == KnotVector{Int}([1,3,2])::KnotVector{Int}
        @test k1 != k2
        @test k1.vector !== copy(k1).vector

        @test KnotVector{Int}([1,2]) isa KnotVector{Int}
        @test KnotVector{Int}([1,2]) isa KnotVector{Int}
        @test KnotVector{Int}([1,2.]) isa KnotVector{Int}
        @test KnotVector{Int}(k1) isa KnotVector{Int}
        @test KnotVector(k1) isa KnotVector{Int}
        @test AbstractKnotVector{Float64}(k1) isa KnotVector{Float64}

        @test knotvector"11111" == KnotVector([1, 2, 3, 4, 5])
        @test knotvector"123" == KnotVector([1, 2, 2, 3, 3, 3])
        @test knotvector" 2 2 2" == KnotVector([2, 2, 4, 4, 6, 6])
        @test knotvector"020202" == KnotVector([2, 2, 4, 4, 6, 6])
        @test knotvector"     1" == KnotVector([6])
    end

    @testset "eltype" begin
        @test eltype(KnotVector([1,2,3])) == Int
        @test eltype(KnotVector{Int}([1,2,3])) == Int
        @test eltype(KnotVector{Real}([1,2,3])) == Real
        @test eltype(KnotVector{Float64}([1,2,3])) == Float64
        @test eltype(KnotVector{BigInt}([1,2,3])) == BigInt
        @test eltype(KnotVector{Rational{Int}}([1,2,3])) == Rational{Int}
    end

    @testset "_vec" begin
        @test BasicBSpline._vec(KnotVector([1,2,3])) isa Vector{Int}
        @test BasicBSpline._vec(EmptyKnotVector()) isa SVector{0,Bool}
    end

    @testset "EmptyKnotVector" begin
        k = EmptyKnotVector()
        @test k === copy(k)
    end

    @testset "zeros" begin
        @test KnotVector(Float64[]) == zero(KnotVector)
        @test KnotVector(Float64[]) == 0*k1 == k1*0 == zero(k1)
        @test KnotVector(Float64[]) == EmptyKnotVector()
        @test KnotVector(Float64[]) |> isempty
        @test KnotVector(Float64[]) |> iszero
        # @test KnotVector() this shold be return error
        @test EmptyKnotVector() |> isempty
        @test EmptyKnotVector() |> iszero
        @test EmptyKnotVector() == KnotVector(Float64[])
        @test KnotVector(Float64[]) == EmptyKnotVector()
        @test EmptyKnotVector{Bool}() === EmptyKnotVector() === zero(EmptyKnotVector) === zero(EmptyKnotVector{Bool})
        @test EmptyKnotVector{Int}() == EmptyKnotVector()
        @test EmptyKnotVector{Int}() == EmptyKnotVector{Int}()
        @test EmptyKnotVector{Int}() == EmptyKnotVector{Float64}()
        @test EmptyKnotVector() == EmptyKnotVector{Int}()
        @test EmptyKnotVector{Int}() == EmptyKnotVector{Int}()
        @test EmptyKnotVector{Float64}() == EmptyKnotVector{Int}()
        @test EmptyKnotVector{Int}() === zero(EmptyKnotVector{Int}())
        @test EmptyKnotVector{Int}() === zero(AbstractKnotVector{Int})
        @test EmptyKnotVector{Int}() !== zero(AbstractKnotVector)
        @test EmptyKnotVector{Bool}() === zero(AbstractKnotVector)
    end

    @testset "length" begin
        @test length(KnotVector(Float64[])) == 0
        @test length(KnotVector([1,2,2,3])) == 4
    end

    @testset "iterator" begin
        for t in k1
            @test t in k2
        end
        @test k1[1] == 1
        @test k2[end] == 3
        @test collect(k2) isa Vector{Int}
        @test [k2...] isa Vector{Int}
        @test collect(k2) == [k2...]
        @test collect(k2) != k2
    end

    @testset "addition, multiply" begin
        @test KnotVector([-1,2,3]) + 2 * KnotVector([2,5]) == KnotVector([-1,2,2,2,3,5,5])
        @test KnotVector([-1,2,3]) + KnotVector([2,5]) * 2 == KnotVector([-1,2,2,2,3,5,5])
        @test k1 + k3 == KnotVector([1,2,2,3,4,5])
        @test 2 * KnotVector([2,3]) == KnotVector([2,2,3,3])

        # EmptyKnotVector
        _k1 = k1 + EmptyKnotVector()
        _k2 = k2 + EmptyKnotVector{Int}()
        _k3 = k3 + EmptyKnotVector{Float64}()
        @test _k1.vector === k1.vector
        @test _k2.vector === k2.vector
        @test _k3.vector !== k3.vector
        _k1 = EmptyKnotVector() + k1
        _k2 = EmptyKnotVector{Int}() + k2
        _k3 = EmptyKnotVector{Float64}() + k3
        @test _k1.vector === k1.vector
        @test _k2.vector === k2.vector
        @test _k3.vector !== k3.vector
        @test EmptyKnotVector{Int}() + EmptyKnotVector{Bool}() isa EmptyKnotVector{Int}
        @test EmptyKnotVector{Int}()*0 === EmptyKnotVector{Int}()
        @test EmptyKnotVector{Float64}()*1 === EmptyKnotVector{Float64}()
        @test EmptyKnotVector{BigFloat}()*2 === EmptyKnotVector{BigFloat}()
        @test_throws DomainError EmptyKnotVector()*(-1)
        @test (EmptyKnotVector()
                === EmptyKnotVector() + EmptyKnotVector()
                === EmptyKnotVector() * 2
                === 2 * EmptyKnotVector())
        @test EmptyKnotVector() isa EmptyKnotVector{Bool}

        # type promotion
        @test KnotVector{Int}([1,2]) + KnotVector([3]) == KnotVector([1,2,3])
        @test KnotVector{Int}([1,2]) + KnotVector([3]) isa KnotVector{Int}
        @test KnotVector{Int}([1,2]) + KnotVector{Rational{Int}}([3]) == KnotVector([1,2,3])
        @test KnotVector{Int}([1,2]) + KnotVector{Rational{Int}}([3]) isa KnotVector{Rational{Int}}
        @test KnotVector{Int}([1,2]) * 0 == KnotVector(Float64[])
        @test KnotVector{Int}([1,2]) * 0 isa KnotVector{Int}
        @test KnotVector{Int}(Int[]) isa KnotVector{Int}
        @test EmptyKnotVector{Irrational{:π}}() + EmptyKnotVector{Rational{Int}}() isa EmptyKnotVector{Float64}
    end

    @testset "unique" begin
        @test unique(k1) == k1
        @test unique(k2) == k1
        @test unique(EmptyKnotVector()) === EmptyKnotVector()
    end

    @testset "inclusive relation" begin
        k1 = KnotVector([1,2,3])
        k2 = KnotVector([1,2,2,3])
        k3 = KnotVector([1,2,2,3,5])
        k4 = KnotVector([1,2,3,5])
        k5 = KnotVector([1,2,3,4])
        k6 = EmptyKnotVector()
        k7 = EmptyKnotVector{Real}()
        k8 = EmptyKnotVector{Float64}()

        @test k1 ⊆ k1
        @test k1 ⊇ k1
        @test k2 ⊆ k3
        @test k2 ⊈ k4
        @test k3 ⊇ k2
        @test k4 ⊉ k2

        @test !(k5 ⊊ k1)
        @test !(k1 ⊊ k1)
        @test k1 ⊊ k5
        @test !(k1 ⊋ k5)
        @test !(k1 ⊋ k1)
        @test k5 ⊋ k1

        @test k6 ⊆ k7 ⊆ k8 ⊆ k1
        @test k6 ⊆ k7 ⊆ k8 ⊆ k2
        @test k6 ⊆ k7 ⊆ k8 ⊆ k3
        @test k6 ⊆ k7 ⊆ k8 ⊆ k4
        @test k6 ⊆ k7 ⊆ k8 ⊆ k5
    end

    @testset "string" begin
        k = KnotVector([1,2,2,3])
        @test string(k) == "KnotVector([1, 2, 2, 3])"
        @test string(KnotVector(Float64[])) == "KnotVector(Float64[])"
        @test string(KnotVector(Int[])) == "KnotVector(Int64[])"

        k1 = UniformKnotVector(1:1:3)
        k2 = UniformKnotVector(Base.OneTo(3))
        k3 = UniformKnotVector(1:4)
        @test string(k1) == "UniformKnotVector(1:1:3)"
        @test string(k2) == "UniformKnotVector(Base.OneTo(3))"
        @test string(k3) == "UniformKnotVector(1:4)"

        @test string(EmptyKnotVector()) == "EmptyKnotVector{Bool}()"
        @test string(EmptyKnotVector{Int}()) == "EmptyKnotVector{Int64}()"
    end

    @testset "float" begin
        k = KnotVector([1,2,3])
        @test float(k) isa KnotVector{Float64}
        @test k == float(k)

        k = KnotVector{Float32}([1,2,3])
        @test float(k) isa KnotVector{Float32}
        @test k == float(k)

        k = KnotVector{Rational{Int}}([1,2,3])
        @test float(k) isa KnotVector{Float64}
        @test k == float(k)
    end

    @testset "other operators" begin
        k = KnotVector([1,2,2,3])
        @test countknots(k, 0.3) == 0
        @test countknots(k, 1.0) == 1
        @test countknots(k, 2.0) == 2
        @test 1 ∈ k
        @test 1.5 ∉ k

        k = KnotVector([-2.0, -1.0, -0.0, -0.0, 0.0, 1.0, 3.0])
        @test countknots(k,0) == 3
        @test countknots(k,0.0) == 3
        @test countknots(k,-0.0) == 3
    end
end
