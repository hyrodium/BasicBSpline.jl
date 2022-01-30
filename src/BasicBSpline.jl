module BasicBSpline

using LinearAlgebra
using IntervalSets
using StaticArrays
using FastGaussQuadrature

export AbstractKnotVector, KnotVector, UniformKnotVector
export BSplineSpace, UniformBSplineSpace
export dim, ⊑, ⊒, ⋢, ⋣, ≃
export BSplineDerivativeSpace
export BSplineManifold, CustomBSplineManifold
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, domain, changebasis, expandspace
export bsplinebasisall, intervalindex
export refinement, bsplinespaces, controlpoints
export isdegenerate, isnondegenerate, exactdim
export degree, knotvector
export 𝔫
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
include("_CustomBSplineManifold.jl")
include("_Refinement.jl")
include("_Fitting.jl")

end # module
