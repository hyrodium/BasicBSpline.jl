using BasicBSpline
using IntervalSets
using LinearAlgebra
using Test
using Random

@testset "BasicBSpline.jl" begin
    include("test_Knots.jl")
    include("test_BSplineSpace.jl")
    include("test_BSplineBasis.jl")
    include("test_ChangeBasis.jl")
    include("test_FastBSplineSpace.jl")
    include("test_BSplineManifold.jl")
    include("test_Integral.jl")
    include("test_Fitting.jl")
end
