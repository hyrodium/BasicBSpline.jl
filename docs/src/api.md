# API

# Public
```@docs
KnotVector
UniformKnotVector
SubKnotVector
EmptyKnotVector
AbstractKnotVector
length(k::AbstractKnotVector)
Base.:+(k1::KnotVector{T}, k2::KnotVector{T}) where T
*(m::Integer, k::AbstractKnotVector)
Base.issubset(k::KnotVector, k′::KnotVector)
unique(k::AbstractKnotVector)
countknots(k::AbstractKnotVector, t::Real)
@knotvector_str
bsplinebasis₊₀
bsplinebasis₋₀
bsplinebasis
BasicBSpline.bsplinebasis₋₀I
intervalindex
bsplinebasisall
bsplinesupport
dim
exactdim_R(P::BSplineSpace)
exactdim_I(P::BSplineSpace)
BSplineSpace
isnondegenerate
isdegenerate(P::BSplineSpace)
BSplineManifold
unbounded_mapping
RationalBSplineManifold
```

```@docs
BasicBSpline.refinement_R
BasicBSpline.refinement_I
```

```@docs
BasicBSplineFitting.innerproduct_R
BasicBSplineFitting.innerproduct_I
```

## Private
Note that the following methods are considered private methods, and changes in their behavior are not considered breaking changes.

```@docs
BasicBSpline.r_nomial
BasicBSpline._vec
BasicBSpline._lower_R
BasicBSpline._changebasis_R
BasicBSpline._changebasis_I
BasicBSpline.__changebasis_I
BasicBSpline._changebasis_sim
```
