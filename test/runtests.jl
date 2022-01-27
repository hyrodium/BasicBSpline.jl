using BasicBSpline
using IntervalSets
using LinearAlgebra
using Test
using Random
using StaticArrays
using GeometryBasics

@testset "BasicBSpline.jl" begin
    include("test_KnotVector.jl")
    include("test_BSplineSpace.jl")
    include("test_BSplineBasis.jl")
    include("test_UniformKnotVector.jl")
    include("test_UniformBSplineSpace.jl")
    include("test_Derivative.jl")
    include("test_ChangeBasis.jl")
    include("test_BSplineManifold.jl")
    include("test_Fitting.jl")
end
