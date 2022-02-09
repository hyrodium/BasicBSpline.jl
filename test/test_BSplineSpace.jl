@testset "BSplineSpace" begin
    Random.seed!(42)

    @testset "dimension, dengenerate" begin
        P1 = BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9]))
        @test bsplinesupport(P1, 2) == 3..8
        @test dim(P1) == 3
        @test exactdim(P1) == 3
        @test isnondegenerate(P1) == true
        @test isdegenerate(P1) == false

        P2 = BSplineSpace{1}(KnotVector([1, 3, 3, 3, 8, 9]))
        @test isnondegenerate(P2) == false
        @test isdegenerate(P2) == true
        @test isdegenerate(P2,1) == false
        @test isdegenerate(P2,2) == true
        @test isdegenerate(P2,3) == false
        @test isdegenerate(P2,4) == false
        @test isnondegenerate(P2,1) == true
        @test isnondegenerate(P2,2) == false
        @test isnondegenerate(P2,3) == true
        @test isnondegenerate(P2,4) == true
        @test dim(P2) == 4
        @test exactdim(P2) == 3
    end

    @testset "denegeration with lower degree" begin
        k = KnotVector([1, 3, 3, 3, 6, 8, 9])
        @test isnondegenerate(BSplineSpace{2}(k))
        @test isdegenerate(BSplineSpace{1}(k))
    end

    @testset "subset" begin
        P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
        P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
        P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
        @test P1 ⊆ P2
        @test P1 ⋢ P2
        @test P1 ⊆ P3
        @test P2 ⊈ P3
    end

    @testset "equality" begin
        P1 = BSplineSpace{2}(KnotVector([1, 2, 3, 5, 8, 8, 9]))
        P2 = BSplineSpace{2}(KnotVector([1, 2, 3, 5, 8, 8, 9]))
        P3 = BSplineSpace{1}(KnotVector([1, 2, 3, 5, 8, 8, 9]))
        P4 = BSplineSpace{1}(KnotVector([1, 2, 3, 5, 8, 8, 9]))

        @test P1 ⊆ P2
        @test P2 ⊆ P1
        @test P1 == P2

        @test P3 ⊆ P4
        @test P4 ⊆ P3
        @test P3 == P4

        @test P1 != P3
    end

    @testset "sqsubset and lowered B-spline space" begin
        P4 = BSplineSpace{1}(KnotVector([1, 2, 3, 4, 5]))
        P5 = BSplineSpace{2}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
        _P4 = BSplineSpace{degree(P4)-1}(knotvector(P4)[2:end-1])
        _P5 = BSplineSpace{degree(P5)-1}(knotvector(P5)[2:end-1])
        @test P4 ⊑ P4
        @test P4 ⊑ P5
        @test P5 ⊒ P4
        @test P5 ⋢ P4
        @test P4 ⋣ P5
        @test (P4 ⊑ P4) == (_P4 ⊑ _P4)
        @test (P4 ⊑ P5) == (_P4 ⊑ _P5)
        @test (P5 ⊑ P4) == (_P5 ⊑ _P4)
    end

    @testset "subset but not sqsubset" begin
        P6 = BSplineSpace{2}(KnotVector(1:8))
        P7 = BSplineSpace{2}(KnotVector(1:10))
        @test P6 ⊆ P6
        @test P7 ⊆ P7
        @test P6 ⊆ P7
        @test P7 ⊈ P6

        @test P6 ⊆ P6
        @test P7 ⊆ P7
        @test P6 ⊆ P7
        @test P7 ⊈ P6
        @test domain(P6) ⊆ domain(P7)
    end

    @testset "sqsubset but not subset" begin
        P6 = BSplineSpace{2}(KnotVector(1,2,3,4,5,6,7))
        P7 = BSplineSpace{2}(KnotVector(0,1,3,4,5,8,9))
        @test P6 ⊑ P7
        @test P7 ⊑ P6
        @test P6 ≃ P7
        @test P6 ⊈ P7
        @test P7 ⊈ P6
        @test domain(P6) == domain(P7)
    end
end
