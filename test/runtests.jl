using BasicBSpline
using ChainRulesTestUtils
using ChainRulesCore
using IntervalSets
using LinearAlgebra
using SparseArrays
using Test
using Random
using StaticArrays
using GeometryBasics
using Plots
using Aqua
using BasicBSpline: clampknotvector, clampknotvector!, isclamped

if VERSION ≥ v"1.9.0"
    # Disable ambiguities tests for ChainRulesCore.frule
    # TODO enable unbound_args again
    Aqua.test_all(BasicBSpline; ambiguities=false, unbound_args=false)
    Aqua.test_ambiguities(BasicBSpline; exclude=[ChainRulesCore.frule])
end
Random.seed!(42)

include("test_util.jl")
include("test_KnotVector.jl")
include("test_BSplineSpace.jl")
include("test_BSplineBasis.jl")
include("test_UniformKnotVector.jl")
include("test_UniformBSplineSpace.jl")
include("test_UniformBSplineBasis.jl")
include("test_Derivative.jl")
include("test_ChangeBasis.jl")
include("test_BSplineManifold.jl")
include("test_RationalBSplineManifold.jl")
include("test_Refinement.jl")
include("test_ChainRules.jl")
include("test_Plots.jl")
