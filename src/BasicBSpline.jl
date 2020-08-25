module BasicBSpline

using IntervalSets
using EllipsisNotation
using FastGaussQuadrature
using Statistics
using LinearAlgebra

export Knots, ♯, BSplineSpace, 𝒫, dim, ⊑
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, bsplineunity, changebasis
export BSplineManifold, refinement, mapping
export isproper, properdim
export degree, knots
export 𝔫
export FastBSplineSpace, FastBSplineManifold
export AbstractBSplineManifold, AbstractBSplineSpace
export FittingControlPoints

const MAX_DEGREE = 3

include("_Knots.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_FastBSplineSpace.jl")
include("_FastBSplineBasis.jl")
include("_BSplineManifold.jl")
include("_Refinement.jl")
include("_Fitting.jl")

end # module
