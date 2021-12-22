module BasicBSpline

using LinearAlgebra
using IntervalSets
using FastGaussQuadrature

export Knots, ♯, BSplineSpace, dim, ⊑, ⊒, ⋢, ⋣, ≃
export BSplineManifold, CustomBSplineManifold
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, bsplineunity, changebasis, lower
export bsplinebasisall, intervalindex
export refinement, bsplinespaces, controlpoints
export isproper, properdim
export degree, knots
export 𝔫
export BSplineCurve, BSplineSurface, BSplineSolid
export AbstractBSplineManifold, AbstractBSplineSpace
export fittingcontrolpoints

const MAX_DEGREE = 5

include("_Knots.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_ChangeBasis.jl")
include("_BSplineManifold.jl")
include("_CustomBSplineManifold.jl")
include("_Refinement.jl")
include("_Integral.jl")
include("_Fitting.jl")

end # module
