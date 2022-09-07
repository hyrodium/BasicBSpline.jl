using BasicBSpline
import BasicBSpline.EmptyKnotVector
using ChainRulesTestUtils
using IntervalSets
using LinearAlgebra
using Test
using Random
using StaticArrays
using GeometryBasics
using Plots
using Aqua

Aqua.test_all(BasicBSpline)

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
include("test_Fitting.jl")
include("test_ChainRules.jl")
include("test_Plots.jl")
