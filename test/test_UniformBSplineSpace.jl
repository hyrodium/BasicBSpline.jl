@testset "UniformBSplineSpace" begin
    Random.seed!(42)

    k1 = UniformKnotVector(1:8)
    P1 = UniformBSplineSpace{2}(k1)
    P1′ = BSplineSpace{2}(k1)
    @test bsplinesupport(P1, 2) == bsplinesupport(P1′, 2) == 2..5
    @test dim(P1) == dim(P1′) == 5
    @test exactdim(P1) == exactdim(P1′) == 5
    @test isnondegenerate(P1) == true
    @test isdegenerate(P1) == false
    @test P1 == P1′

    k2 = UniformKnotVector(range(1,1,length=8))
    P2 = UniformBSplineSpace{2}(k2)
    P2′ = BSplineSpace{2}(k2)
    @test bsplinesupport(P2, 2) == bsplinesupport(P2′, 2) == 1..1
    @test dim(P2) == dim(P2′) == 5
    @test exactdim(P2) == exactdim(P2′) == 0
    @test isnondegenerate(P2) == false
    @test isdegenerate(P2) == true
    @test P2 == P2′

    k3 = UniformKnotVector(1.0:0.5:12.0)
    P3 = UniformBSplineSpace{2}(k3)
    @test P1 ⊆ P3
    @test P2 ⋢ P3
    @test P1 ⋢ P3
    @test P2 ⋢ P3

    k4 = UniformKnotVector(1.0:0.5:12.0)
    P3 = UniformBSplineSpace{2}(k3)
    @test P1 ⊆ P3
    @test P2 ⋢ P3
    @test P1 ⋢ P3
    @test P2 ⋢ P3

    # P2 = BSplineSpace{1}(KnotVector([1, 3, 3, 3, 8, 9]))
    # @test isnondegenerate(P2) == false
    # @test isdegenerate(P2) == true
    # @test isdegenerate(P2,1) == false
    # @test isdegenerate(P2,2) == true
    # @test isdegenerate(P2,3) == false
    # @test isdegenerate(P2,4) == false
    # @test isnondegenerate(P2,1) == true
    # @test isnondegenerate(P2,2) == false
    # @test isnondegenerate(P2,3) == true
    # @test isnondegenerate(P2,4) == true
    # @test dim(P2) == 4
    # @test exactdim(P2) == 3

    # i = 2
    # k = KnotVector([5, 12, 13, 13, 14])
    # p = 2
    # P = BSplineSpace{p}(k)
    # @test bsplinesupport(P) == [5..13, 12..14]
    # @test bsplinesupport(P, i) == 12..14

    # @test isnondegenerate(BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9])))
    # @test isdegenerate(BSplineSpace{1}(KnotVector([1, 3, 3, 3, 8, 9])))

    # @test dim(BSplineSpace{2}(KnotVector([1, 3, 5, 6, 8, 9]))) == 3

    # P1 = BSplineSpace{1}(KnotVector([1, 3, 5, 8]))
    # P2 = BSplineSpace{1}(KnotVector([1, 3, 5, 6, 8, 9]))
    # P3 = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
    # @test P1 ⊆ P2
    # @test P1 ⋢ P2
    # @test P1 ⊆ P3
    # @test P2 ⊈ P3

    # P4 = BSplineSpace{1}(KnotVector([1, 2, 3, 4, 5]))
    # P5 = BSplineSpace{2}(KnotVector([-1, 0.3, 2, 3, 3, 4, 5.2, 6]))
    # _P4 = BSplineSpace{degree(P4)-1}(knotvector(P4)[2:end-1])
    # _P5 = BSplineSpace{degree(P5)-1}(knotvector(P5)[2:end-1])
    # @test P4 ⊑ P4
    # @test P4 ⊑ P5
    # @test P5 ⊒ P4
    # @test P5 ⋢ P4
    # @test P4 ⋣ P5
    # @test (P4 ⊑ P4) == (_P4 ⊑ _P4)
    # @test (P4 ⊑ P5) == (_P4 ⊑ _P5)
    # @test (P5 ⊑ P4) == (_P5 ⊑ _P4)

    # P6 = BSplineSpace{2}(KnotVector(1,2,3,4,5,6,7))
    # P7 = BSplineSpace{2}(KnotVector(0,1,3,4,5,8,9))
    # @test P6 ⊑ P7
    # @test P7 ⊑ P6
    # @test P6 ≃ P7
    # @test P6 ⊈ P7
    # @test P7 ⊈ P6
end
