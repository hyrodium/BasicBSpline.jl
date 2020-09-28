module BasicBSpline

using LinearAlgebra
using EllipsisNotation
using IntervalSets
using FastGaussQuadrature

export Knots, ♯, BSplineSpace, 𝒫, dim, ⊑, ⊒, ⋢, ⋣, ≃
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, bsplineunity, changebasis
export BSplineManifold, refinement, mapping, bsplinespaces, controlpoints
export isproper, properdim
export degree, knots
export 𝔫
export FastBSplineSpace, FastBSplineManifold
export BSplineCurve, BSplineSurface, BSplineSolid
export AbstractBSplineManifold, AbstractBSplineSpace
export fittingcontrolpoints

const MAX_DEGREE = 3

include("_Knots.jl")
include("_BSplineSpace.jl")
include("_FastBSplineSpace.jl")
include("_BSplineBasis.jl")
include("_FastBSplineBasis.jl")
include("_BSplineManifold.jl")
include("_FastBSplineManifold.jl")
include("_Refinement.jl")
include("_Integral.jl")
include("_Fitting.jl")

end # module
