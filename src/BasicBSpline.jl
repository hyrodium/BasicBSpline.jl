module BasicBSpline

using LinearAlgebra
using IntervalSets
using FastGaussQuadrature

export KnotVector, BSplineSpace, dim, ⊑, ⊒, ⋢, ⋣, ≃
export BSplineDerivativeSpace
export BSplineManifold, CustomBSplineManifold
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, domain, changebasis, expandspace
export bsplinebasisall, intervalindex
export refinement, bsplinespaces, controlpoints
export isproper, properdim
export degree, knots
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
