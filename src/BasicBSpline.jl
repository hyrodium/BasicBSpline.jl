module BasicBSpline

using LinearAlgebra
using IntervalSets
using StaticArrays
using FastGaussQuadrature

# Types
export AbstractKnotVector, KnotVector, UniformKnotVector
export AbstractBSplineSpace, BSplineSpace, UniformBSplineSpace
export AbstractBSplineManifold, BSplineManifold, RationalBSplineManifold
export BSplineDerivativeSpace

# B-spline basis functions
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinebasisall, intervalindex

# B-spline space
export dim, exactdim, ⊑, ⊒, ⋢, ⋣, ≃
export domain, changebasis
export bsplinesupport, bsplinesupport_R, bsplinesupport_I
export expandspace, expandspace_R, expandspace_I
export isdegenerate, isdegenerate_R, isdegenerate_I
export isnondegenerate, isnondegenerate_R, isnondegenerate_I
export countknots
export degree

# Getter methods
export bsplinespaces, controlpoints, weights
export knotvector
export bsplinespace

# Useful functions
export refinement
export fittingcontrolpoints

include("_util.jl")
include("_KnotVector.jl")
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

end # module
