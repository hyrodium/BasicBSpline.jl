module BasicBSpline

using LinearAlgebra
using IntervalSets
using StaticArrays
using FastGaussQuadrature

export KnotVector, AbstractKnotVector
export BSplineSpace, dim, ⊑, ⊒, ⋢, ⋣, ≃
export BSplineDerivativeSpace
export BSplineManifold, CustomBSplineManifold
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, domain, changebasis, expandspace
export bsplinebasisall, intervalindex
export refinement, bsplinespaces, controlpoints
export isdegenerate, isnondegenerate, properdim
export degree, knotvector
export 𝔫
export AbstractBSplineManifold, AbstractBSplineSpace
export fittingcontrolpoints

include("_KnotVector.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_DerivativeSpace.jl")
include("_DerivativeBasis.jl")
include("_ChangeBasis.jl")
include("_BSplineManifold.jl")
include("_CustomBSplineManifold.jl")
include("_Refinement.jl")
include("_Integral.jl")
include("_Fitting.jl")

end # module
