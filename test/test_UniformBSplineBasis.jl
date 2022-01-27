@testset "UniformKnotVector" begin
    @testset "bsplinebasis" begin
        k = KnotVector(0:10)
        k1 = UniformKnotVector(0:10)
        k2 = UniformKnotVector(0:0.5:5)
        k3 = UniformKnotVector(0:2:20)
        @test eltype(k1) == Int
        for p in 0:3
            P = BSplineSpace{p}(k)
            P1 = UniformBSplineSpace{p}(k1)
            P2 = UniformBSplineSpace{p}(k2)
            P3 = UniformBSplineSpace{p}(k3)
            @test dim(P) == dim(P1)
            @test dim(P) == dim(P2)
            @test dim(P) == dim(P3)
            n = dim(P)
            for _ in 1:20
                t = 10*rand()
                for i in 1:n
                    @test bsplinebasis(P1,i,t) ≈ bsplinebasis(P,i,t)
                    @test bsplinebasis(P2,i,t/2) ≈ bsplinebasis(P,i,t)
                    @test bsplinebasis(P3,i,t*2) ≈ bsplinebasis(P,i,t)
                    @test (t ∈ bsplinesupport(P1,i)) == (bsplinebasis(P1,i,t) != 0)
                    @test (t/2 ∈ bsplinesupport(P2,i)) == (bsplinebasis(P2,i,t/2) != 0)
                    @test (t*2 ∈ bsplinesupport(P3,i)) == (bsplinebasis(P3,i,t*2) != 0)
                    @test bsplinebasis(P1,i,t) isa Float64
                    @test bsplinebasis(P2,i,t) isa Float64
                    @test bsplinebasis(P3,i,t) isa Float64
                end
            end
        end
    end
end
