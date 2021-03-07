@testset "Knots" begin
    @test Knots([1,2,3]) == Knots(1:3) == Knots(1,2,3)
    @test Knots(2,4,5) == Knots([2,4,5])
    @test zero(Knots) == Knots() == Knots([])
    @test Knots(1:3) == Knots([3,2,1])

    k = Knots([1,2,2,3])
    @test ♯(k) == length(k) == 4

    @test Knots([-1,2,3]) + 2 * Knots([2,5]) == Knots([-1,2,2,2,3,5,5])
    @test Knots([1,2,3]) + Knots([2,4,5]) == Knots([1,2,2,3,4,5])
    @test 2 * Knots([2,3]) == Knots([2,2,3,3])

    @test unique(Knots([1,2,2,3])) == Knots([1,2,3])

    @test Knots([1,2,3]) ⊆ Knots([1,2,3])
    @test Knots([1,2,3]) ⊇ Knots([1,2,3])
    @test Knots([1,2,2,3]) ⊆ Knots([1,2,2,3,5])
    @test Knots([1,2,2,3]) ⊈ Knots([1,2,3,5])
    @test Knots([1,2,2,3,5]) ⊇ Knots([1,2,2,3])
    @test Knots([1,2,3,5]) ⊉ Knots([1,2,2,3])

    @test 𝔫(k, 0.3) == 0
    @test 𝔫(k, 1.0) == 1
    @test 𝔫(k, 2.0) == 2
    @test 1 ∈ k
    @test 1.5 ∉ k
    @test BasicBSpline._knotindex(k,1.5) == 1
    @test BasicBSpline._knotindex(k,2.5) == 3
    @test BasicBSpline._knotindex₋₀(k,1.5) == 1
    @test BasicBSpline._knotindex₋₀(k,2.5) == 3
    @test BasicBSpline._knotindex₊₀(k,1.5) == 1
    @test BasicBSpline._knotindex₊₀(k,2.5) == 3

    @test BasicBSpline._knotindex(k,1) == 1
    @test BasicBSpline._knotindex(k,2) == 3
    @test BasicBSpline._knotindex(k,3) == 3
    @test BasicBSpline._knotindex₋₀(k,1) == nothing
    @test BasicBSpline._knotindex₋₀(k,2) == 1
    @test BasicBSpline._knotindex₋₀(k,3) == 3
    @test BasicBSpline._knotindex₊₀(k,1) == 1
    @test BasicBSpline._knotindex₊₀(k,2) == 3
    @test BasicBSpline._knotindex₊₀(k,3) == nothing

    @test string(k) == "Knots([1.0, 2.0, 2.0, 3.0])"
    @test string(Knots()) == "Knots([])"
end
