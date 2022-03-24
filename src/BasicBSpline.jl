module BasicBSpline

using LinearAlgebra
using IntervalSets
using StaticArrays
using FastGaussQuadrature

export AbstractKnotVector, KnotVector, UniformKnotVector
export BSplineSpace, UniformBSplineSpace
export dim, ⊑, ⊒, ⋢, ⋣, ≃
export BSplineDerivativeSpace
export BSplineManifold, RationalBSplineManifold
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export domain, changebasis
export bsplinesupport, bsplinesupport_R, bsplinesupport_I
export expandspace, expandspace_R, expandspace_I
export isdegenerate, isdegenerate_R, isdegenerate_I
export isnondegenerate, isnondegenerate_R, isnondegenerate_I
export bsplinebasisall, intervalindex
export bsplinespaces, controlpoints, weights
export refinement
export exactdim
export degree, knotvector
export countknots
export AbstractBSplineManifold, AbstractBSplineSpace
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
