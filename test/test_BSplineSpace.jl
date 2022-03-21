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

        P3 = BSplineSpace{2}(KnotVector(0,0,0,1,1,1,2))
        @test dim(P3) == 4
        @test exactdim(P3) == 4
        @test isnondegenerate(P3) == true
        @test isdegenerate(P3) == false
        @test isdegenerate(P3,1) == false
        @test isdegenerate(P3,2) == false
        @test isdegenerate(P3,3) == false
        @test isdegenerate(P3,4) == false
        @test isnondegenerate(P3,1) == true
        @test isnondegenerate(P3,2) == true
        @test isnondegenerate(P3,3) == true
        @test isnondegenerate(P3,4) == true
        @test isnondegenerate_I(P3) == false
        @test isdegenerate_I(P3) == true
        @test isdegenerate_I(P3,1) == true
        @test isdegenerate_I(P3,2) == true
        @test isdegenerate_I(P3,3) == true
        @test isdegenerate_I(P3,4) == false
        @test isnondegenerate_I(P3,1) == false
        @test isnondegenerate_I(P3,2) == false
        @test isnondegenerate_I(P3,3) == false
        @test isnondegenerate_I(P3,4) == true
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

    @testset "expandspace" begin
        P1 = BSplineSpace{1}(KnotVector(0,0,0,0,0))
        P2 = BSplineSpace{3}(KnotVector(1:8))
        P3 = BSplineSpace{2}(KnotVector(1:8)+3*KnotVector(0))
        P4 = BSplineSpace{2}(KnotVector(0,0,0,1,1,1,2))

        @test P1 == expandspace_R(P1) == expandspace_I(P1) == expandspace(P1)
        @test P2 == expandspace_R(P2) == expandspace_I(P2) == expandspace(P2)
        @test P3 == expandspace_R(P3) == expandspace_I(P3) == expandspace(P3)
        @test P4 == expandspace_R(P4) == expandspace_I(P4) == expandspace(P4)

        @test P1 ⊆ expandspace_R(P1, p₊=1)
        @test P2 ⊆ expandspace_R(P2, p₊=1)
        @test P3 ⊆ expandspace_R(P3, p₊=1)
        @test P4 ⊆ expandspace_R(P4, p₊=1)

        @test_broken P1 ⊑ expandspace_I(P1, p₊=1) == expandspace(P1, p₊=1)
        @test P2 ⊑ expandspace_I(P2, p₊=1) == expandspace(P2, p₊=1)
        @test P3 ⊑ expandspace_I(P3, p₊=1) == expandspace(P3, p₊=1)
        @test_broken P4 ⊑ expandspace_I(P4, p₊=1) == expandspace(P4, p₊=1)
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
