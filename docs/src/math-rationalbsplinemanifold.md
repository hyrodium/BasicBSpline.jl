# Rational B-spline manifold

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

[Non-uniform rational basis spline (NURBS)](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline) is also supported in BasicBSpline.jl package.

## Rational B-spline manifold
Rational B-spline manifold is a parametric representation of a shape.

!!! tip "Def.  Rational B-spline manifold"
    For given ``d``-dimensional B-spline basis functions ``B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)}``, given points ``\bm{a}_{i^1 \dots i^d} \in V`` and real numbers ``w_{i^1 \dots i^d} > 0``, rational B-spline manifold is defined by the following equality:
    ```math
    \bm{p}(t^1,\dots,t^d; \bm{a}_{i^1 \dots i^d}, w_{i^1 \dots i^d})
    =\sum_{i^1,\dots,i^d}
    \frac{(B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)})(t^1,\dots,t^d) w_{i^1 \dots i^d}}
    {\sum\limits_{j^1,\dots,j^d}(B_{(j^1,p^1,k^1)} \otimes \cdots \otimes B_{(j^d,p^d,k^d)})(t^1,\dots,t^d) w_{j^1 \dots j^d}}
    \bm{a}_{i^1 \dots i^d}
    ```
    Where ``\bm{a}_{i^1,\dots,i^d}`` are called **control points**, and ``w_{i^1 \dots i^d}`` are called **weights**.

```@docs
RationalBSplineManifold
```

## Properties
Similar to `BSplineManifold`, `RationalBSplineManifold` supports the following methods and properties.

* currying
* `refinement`
* Affine commutativity
