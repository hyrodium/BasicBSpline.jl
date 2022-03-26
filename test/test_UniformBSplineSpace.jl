@testset "UniformBSplineSpace" begin
    Random.seed!(42)

    @testset "constructor, subset" begin
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
    end

    @testset "dimension, dengenerate" begin
        P1 = UniformBSplineSpace{2}(UniformKnotVector(1:8))
        @test bsplinesupport(P1, 2) == 2..5
        @test dim(P1) == 5
        @test exactdim(P1) == 5
        @test isnondegenerate(P1) == true
        @test isdegenerate(P1) == false

        P2 = UniformBSplineSpace{2}(UniformKnotVector(range(1,1,length=8)))
        @test isnondegenerate(P2) == false
        @test isdegenerate(P2) == true
        @test isdegenerate(P2,1) == true
        @test isdegenerate(P2,2) == true
        @test isdegenerate(P2,3) == true
        @test isdegenerate(P2,4) == true
        @test isnondegenerate(P2,1) == false
        @test isnondegenerate(P2,2) == false
        @test isnondegenerate(P2,3) == false
        @test isnondegenerate(P2,4) == false
        @test dim(P2) == 5
        @test exactdim(P2) == 0
    end

    @testset "support" begin
        k = UniformKnotVector(0:9)
        P = UniformBSplineSpace{2}(k)
        @test bsplinesupport(P,1) == 0..3
        @test bsplinesupport(P,2) == 1..4
    end

    @testset "equality" begin
        P1 = UniformBSplineSpace{2}(UniformKnotVector(Base.OneTo(8)))
        P2 = UniformBSplineSpace{2}(UniformKnotVector(1:8))
        P3 = UniformBSplineSpace{2}(UniformKnotVector(1:1:8))
        P4 = UniformBSplineSpace{2}(UniformKnotVector(1.0:1.0:8.0))
        P5 = UniformBSplineSpace{1}(UniformKnotVector(1.0:1.0:8.0))
        P6 = UniformBSplineSpace{2}(UniformKnotVector(1.0:2.0:8.0))

        @test P1 == P2 == P3 == P4
        @test P1 != P5
        @test P1 != P6
    end

    @testset "subset but not sqsubset" begin
        P1 = UniformBSplineSpace{2}(UniformKnotVector(0:1:8))
        P2 = UniformBSplineSpace{2}(UniformKnotVector(0.0:0.5:12.0))

        @test P1 ⊆ P1
        @test P2 ⊆ P2
        @test P1 ⊆ P2
        @test P2 ⊈ P1

        @test P1 ⊑ P1
        @test P2 ⊑ P2
        @test P1 ⋢ P2
        @test P2 ⋢ P1
    end

    @testset "sqsubset but not subset" begin
        P1 = UniformBSplineSpace{2}(UniformKnotVector(0:1:8))
        P2 = UniformBSplineSpace{2}(UniformKnotVector(1.0:0.5:7.0))

        @test domain(P1) == domain(P2)
        
        @test P1 ⊆ P1
        @test P2 ⊆ P2
        @test P1 ⊈ P2
        @test P2 ⊈ P1

        @test P1 ⊑ P1
        @test P2 ⊑ P2
        @test P1 ⊑ P2
        @test P2 ⋢ P1
    end
end
