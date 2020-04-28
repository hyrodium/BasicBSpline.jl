module BasicBSpline

using IntervalSets

export Knots, BSplineSpace, 𝒫, dim
export BSplineBasis₊₀, BSplineBasis₋₀, BSplineBasis
export BSplineBasis′₊₀, BSplineBasis′₋₀, BSplineBasis′
export BSplineSupport, BSplineCoefficient
export BSplineManifold, Refinement, Mapping, BSplineSvg
export isproper, properdim
export 𝔫

include("_Knots.jl")
include("_BSplineSpace.jl")
include("_BSplineBasis.jl")
include("_BSplineManifold.jl")
include("_Refinement.jl")


end # module
