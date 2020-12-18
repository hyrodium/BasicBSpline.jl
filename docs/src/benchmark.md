# Benchmark
## Concrete type determination; generality vs speed performance
### B-spline space
There are two subtypes of `AbstractBSplineSpace`.
```julia
julia> subtypes(AbstractBSplineSpace)
2-element Array{Any,1}:
 BSplineSpace
 FastBSplineSpace
```
* `BSplineSpace`
    * General degree of polynomial is supported.
* `FastBSplineSpace`
    * 5 or less degree of polynomial is supported.
    * But much faster than `BSplineSpace`

### B-spline manifold
There are five subtypes of `AbstractBSplineManifold`.
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
    * 5 or less degree of polynomial is supported.
    * General dimensional manifold is supported.
    * Much faster than `BSplineManifold`
* `BSplineCurve`
    * 5 or less degree of polynomial is supported.
    * Only one dimensional manifold (curve) is supported.
    * Much faster than `FastBSplineManifold`
* `BSplineSurface`
    * 5 or less degree of polynomial is supported.
    * Only two dimensional manifold (surface) is supported.
    * Much faster than `FastBSplineManifold`
* `BSplineCurve`
    * 5 or less degree of polynomial is supported.
    * Only three dimensional manifold (solid) is supported.
    * Much faster than `FastBSplineManifold`

## Detailed Benchmarks
### B-spline basis function
In this section, we compare the implementation of [Cox-de Boor algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm).

#### Naive implementation
```julia
using BenchmarkTools

function B(i,p,k,t)
    if p==0
        return float(k[i]≤t<k[i+1])
    else
        return B(i,p-1,k,t)*(t-k[i])/(k[i+p]-k[i]) + B(i+1,p-1,k,t)*(k[i+p+1]-t)/(k[i+p+1]-k[i+1])
    end
end

knotvector = sort(rand(10)) # Vector{Float64} with 10 elements
p = 3
i = 1
t = 0.2

@benchmark B(i,p,knotvector,t)
```
Result: 118.757 ns

#### Use type `BSplineSpace`

```julia
using BenchmarkTools
using BasicBSpline
knotvector = sort(rand(10)) # Vector{Float64} with 10 elements
p = 3
i = 1
t = 0.2
k = Knots(knotvector)
P = BSplineSpace(p,k)
@benchmark bsplinebasis(i,P,t)
```
Result: 223.359 ns

Slower than naive implementation.


#### Use type `FastBSplineSpace`

```julia
using BenchmarkTools
using BasicBSpline
knotvector = sort(rand(10)) # Vector{Float64} with 10 elements
p = 3
i = 1
t = 0.2
k = Knots(knotvector)
P = FastBSplineSpace(p,k)
@benchmark bsplinebasis(i,P,t)
```
Result: 46.722 ns

The fastest.

### B-spline manifold

#### Use type `BSplineManifold`

#### Use type `FastBSplineManifold`

#### Use type `BSplineSurface`
