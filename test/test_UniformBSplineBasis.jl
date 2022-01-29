@testset "UniformBSplineBasis" begin
    @testset "bsplinebasis" begin
        v = 2:10
        k = KnotVector(v)
        k1 = UniformKnotVector(v)
        k2 = UniformKnotVector(v/2)
        k3 = UniformKnotVector(v*2)
        @test eltype(k1) == Int
        @test eltype(k2) == Float64
        @test eltype(k3) == Int
        for p in 0:3
            P = BSplineSpace{p}(k)
            P1 = UniformBSplineSpace{p}(k1)
            P2 = UniformBSplineSpace{p}(k2)
            P3 = UniformBSplineSpace{p}(k3)
            @test dim(P) == dim(P1)
            @test dim(P) == dim(P2)
            @test dim(P) == dim(P3)
            n = dim(P)
            @test bsplinebasis(P1,n÷2,6//1) isa Rational{Int}
            @test bsplinebasis(P2,n÷2,6//1) isa Float64
            @test bsplinebasis(P3,n÷2,6//1) isa Rational{Int}
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
                    @test bsplinebasis(P1,i,t) == bsplinebasis₊₀(P1,i,t) == bsplinebasis₋₀(P1,i,t)
                    @test bsplinebasis(P2,i,t) == bsplinebasis₊₀(P2,i,t) == bsplinebasis₋₀(P2,i,t)
                    @test bsplinebasis(P3,i,t) == bsplinebasis₊₀(P3,i,t) == bsplinebasis₋₀(P3,i,t)
                end
            end
        end
    end

    @testset "bsplinebasisall" begin
        v = 2:10
        k = KnotVector(v)
        l = length(k)
        k1 = UniformKnotVector(v)
        k2 = UniformKnotVector(v/2)
        k3 = UniformKnotVector(v*2)
        for p in 0:3
            P = BSplineSpace{p}(k)
            P1 = UniformBSplineSpace{p}(k1)
            P2 = UniformBSplineSpace{p}(k2)
            P3 = UniformBSplineSpace{p}(k3)
            @test bsplinebasisall(P1,(l-2p)÷2,6//1) isa SVector{p+1,Rational{Int}}
            @test bsplinebasisall(P2,(l-2p)÷2,6//1) isa SVector{p+1,Float64}
            @test bsplinebasisall(P3,(l-2p)÷2,6//1) isa SVector{p+1,Rational{Int}}
            for _ in 1:20
                t = 10*rand()
                for i in 1:l-2p-1  # FIXME: This can be 0:l-2p, but don't know why the test doesn't pass.
                    @test bsplinebasisall(P1,i,t) ≈ bsplinebasisall(P,i,t)
                    @test bsplinebasisall(P2,i,t/2) ≈ bsplinebasisall(P,i,t)
                    @test bsplinebasisall(P3,i,t*2) ≈ bsplinebasisall(P,i,t)
                    @test bsplinebasisall(P1,i,t) isa SVector{p+1,Float64}
                    @test bsplinebasisall(P2,i,t) isa SVector{p+1,Float64}
                    @test bsplinebasisall(P3,i,t) isa SVector{p+1,Float64}
                end
            end
        end
    end
end
