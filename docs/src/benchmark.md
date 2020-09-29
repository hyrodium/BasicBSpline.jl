# Benchmark
## Speed and Type
### B-spline space
There are two subtypes for `AbstractBSplineSpace`.
```julia
julia> subtypes(AbstractBSplineSpace)
2-element Array{Any,1}:
 BSplineSpace
 FastBSplineSpace
```
* `BSplineSpace`
    * General degree of polynomial is supported.
* `FastBSplineSpace`
    * 3 or less degree of polynomial is supported.
    * But much faster than `BSplineSpace`

### B-spline manifold
There are five subtypes for `AbstractBSplineManifold`.
```julia
julia> subtypes(AbstractBSplineManifold)
5-element Array{Any,1}:
 BSplineCurve
 BSplineManifold
 BSplineSolid
 BSplineSurface
 FastBSplineManifold
```
* `BSplineManifold`
    * General degree of polynomial is supported.
    * General dimension of polynomial is supported.
* `FastBSplineManifold`
    * 3 or less degree of polynomial is supported.
    * General dimensional manifold is supported.
    * Much faster than `BSplineManifold`
* `BSplineCurve`
    * 3 or less degree of polynomial is supported.
    * Only one dimensional manifold (curve) is supported.
    * Much faster than `FastBSplineManifold`
* `BSplineSurface`
    * 3 or less degree of polynomial is supported.
    * Only two dimensional manifold (surface) is supported.
    * Much faster than `FastBSplineManifold`
* `BSplineCurve`
    * 3 or less degree of polynomial is supported.
    * Only three dimensional manifold (solid) is supported.
    * Much faster than `FastBSplineManifold`

## Detailed Benchmarks

### Naive implementation
(TODO)

### BasicBSpline
(TODO)
