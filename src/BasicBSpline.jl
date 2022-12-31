module BasicBSpline

using LinearAlgebra
using IntervalSets
using StaticArrays
using RecipesBase
using FastGaussQuadrature
using ChainRulesCore

# Types
export AbstractKnotVector, KnotVector, UniformKnotVector, EmptyKnotVector
export AbstractBSplineSpace, BSplineSpace
export AbstractBSplineManifold, BSplineManifold, RationalBSplineManifold
export BSplineDerivativeSpace

# B-spline basis functions
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinebasisall, intervalindex

# B-spline space
export issqsubset, ⊑, ⊒, ⋢, ⋣, ⋤, ⋥, ≃
export dim, exactdim
export domain, changebasis
export bsplinesupport, bsplinesupport_R, bsplinesupport_I
export expandspace, expandspace_R, expandspace_I
export isdegenerate, isdegenerate_R, isdegenerate_I
export isnondegenerate, isnondegenerate_R, isnondegenerate_I
export countknots
export degree
export unbounded_mapping

# Getter methods
export bsplinespaces, controlpoints, weights
export knotvector
export bsplinespace

# Useful functions
export refinement
export fittingcontrolpoints, fittingcontrolpoints_R, fittingcontrolpoints_I

include("_util.jl")
include("_KnotVector.jl")
include("_EmptyKnotVector.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_UniformKnotVector.jl")
include("_UniformBSplineSpace.jl")
include("_UniformBSplineBasis.jl")
include("_DerivativeSpace.jl")
include("_DerivativeBasis.jl")
include("_ChangeBasis.jl")
include("_BSplineManifold.jl")
include("_RationalBSplineManifold.jl")
include("_Refinement.jl")
include("_Fitting.jl")
include("_ChainRules.jl")
include("_PlotRecipes.jl")

end # module
