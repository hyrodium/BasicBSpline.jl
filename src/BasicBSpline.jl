module BasicBSpline

using LinearAlgebra
using SparseArrays
using IntervalSets
using StaticArrays
using RecipesBase
using ChainRulesCore

# Types
export AbstractKnotVector, KnotVector, UniformKnotVector, EmptyKnotVector, SubKnotVector
export BSplineSpace, BSplineDerivativeSpace, UniformBSplineSpace
export BSplineManifold, RationalBSplineManifold

# B-spline basis functions
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinebasisall, intervalindex

# B-spline space
export issqsubset, ⊑, ⊒, ⋢, ⋣, ⋤, ⋥, ≃
export dim, exactdim_R
export domain
export changebasis, changebasis_R, changebasis_I
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

# Macros
export @knotvector_str

include("_util.jl")
include("_KnotVector.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_DerivativeSpace.jl")
include("_DerivativeBasis.jl")
include("_ChangeBasis.jl")
include("_BSplineManifold.jl")
include("_RationalBSplineManifold.jl")
include("_Refinement.jl")
include("_ChainRules.jl")
include("_PlotRecipes.jl")

end # module
