module BasicBSpline

using IntervalSets
using EllipsisNotation

export Knots, ♯, BSplineSpace, 𝒫, dim
export bsplinebasis₊₀, bsplinebasis₋₀, bsplinebasis
export bsplinebasis′₊₀, bsplinebasis′₋₀, bsplinebasis′
export bsplinesupport, changebasis
export BSplineManifold, refinement, mapping
export isproper, properdim
export degree, knots
export 𝔫
export FastBSplineSpace, FastBSplineManifold

const MAX_DEGREE = 3

include("_Knots.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_FastBSplineSpace.jl")
include("_FastBSplineBasis.jl")
include("_BSplineManifold.jl")
include("_Refinement.jl")

end # module
