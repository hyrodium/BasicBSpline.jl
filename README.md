# BasicBSpline.jl

Basic operations for B-spline functions and related things with julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/dev)
[![Build Status](https://github.com/hyrodium/BasicBSpline.jl/workflows/CI/badge.svg)](https://github.com/hyrodium/BasicBSpline.jl/actions)
[![Coverage](https://codecov.io/gh/hyrodium/BasicBSpline.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hyrodium/BasicBSpline.jl)
[![DOI](https://zenodo.org/badge/258791290.svg)](https://zenodo.org/badge/latestdoi/258791290)

![](docs/src/img/BasicBSplineLogo.png)

## Summary
This package provides basic (mathematical) operations for [B-spline](https://en.wikipedia.org/wiki/B-spline).

* B-spline basis function
* Some operations for knot vector
* B-spline manifold (includes curve, surface and solid)
* Refinement for B-spline manifold
* Fitting control points for B-spline manifold

## Comparison to other julia packages for B-spline
* [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
    * >Currently this package's support is best for B-splines and also supports irregular grids.
    * But seems like no method for B-spline manifold.
* [ApproXD.jl](https://github.com/floswald/ApproXD.jl)
    * Its functions are similar to Interpolations.jl.
* [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
    * Wrapper for the dierckx Fortran library.
    * Only 1-d or 2-d B-spline manifold are supported.
    * 5 or less degree of polynomial are supported.
* **[BasicBSpline.jl](https://github.com/hyrodium/BasicBSpline.jl) (this package)**
    * Any degree of polynomial are supported.
    * Refinement algorithm for B-spline manifold.
    * Fitting algorithm by least squares.

## Installation
Install this package

```julia
] add BasicBSpline
```

To export graphics, use [BasicBSplineExporter.jl](https://github.com/hyrodium/BasicBSplineExporter.jl).

```julia
] add https://github.com/hyrodium/BasicBSplineExporter.jl
```

## Example
### B-spline function

```julia
using BasicBSpline
using Plots

k = KnotVector([0.00,1.50,2.50,5.50,8.00,9.00,9.50,10.0])
P0 = BSplineSpace{0}(k) # 0th degree piecewise polynomial space
P1 = BSplineSpace{1}(k) # 1st degree piecewise polynomial space
P2 = BSplineSpace{2}(k) # 2nd degree piecewise polynomial space
P3 = BSplineSpace{3}(k) # 3rd degree piecewise polynomial space
plot(
    plot([t->bsplinebasis(P0,i,t) for i in 1:dim(P0)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(P1,i,t) for i in 1:dim(P1)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(P2,i,t) for i in 1:dim(P2)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(P3,i,t) for i in 1:dim(P3)], 0, 10, ylims=(0,1), legend=false),
    layout=(2,2),
)
```

![](docs/src/img/cover.png)

Try [interactive graph with Desmos graphing calculator](https://www.desmos.com/calculator/ql6jqgdabs)!

### B-spline manifold
```julia
using BasicBSpline
using BasicBSplineExporter

p = 2 # degree of polynomial
k1 = KnotVector(1:8)     # knot vector
k2 = KnotVector(rand(8)) # knot vector
P1 = BSplineSpace{p}(k1) # B-spline space
P2 = BSplineSpace{p}(k2) # B-spline space
n1 = dim(P1) # dimension of B-spline space
n2 = dim(P2) # dimension of B-spline space
ax = [2*i-6.5+rand() for i in 1:dim(P), j in 1:dim(P)] # x-coordinates of random generated control points
ay = [2*j-6.5+rand() for i in 1:dim(P), j in 1:dim(P)] # y-coordinates of random generated control points
a = [ax;;;ay] # n1 × n2 × 2 dimension array
M = BSplineManifold(a,(P,P)) # Define B-spline manifold
save_png("2dim.png", M) # save image
```
![](docs/src/img/2dim.png)

### Refinement
#### h-refinemnet
```julia
k₊=(KnotVector(3.3,4.2),KnotVector(3.8,3.2,5.3)) # additional knotvectors
M_h = refinement(M,k₊=k₊) # refinement of B-spline manifold
save_png("2dim_h-refinement.png", M_h) # save image
```
![](docs/src/img/2dim_h-refinement.png)

Note that this shape and the last shape are identical.

#### p-refinemnet
```julia
p₊=(1,2) # additional degrees
M_p = refinement(M,p₊=p₊) # refinement of B-spline manifold
save_png("2dim_p-refinement.png", M_p) # save image
```
![](docs/src/img/2dim_p-refinement.png)

Note that this shape and the last shape are identical.

### Fitting B-spline manifold
[Try on Desmos graphing graphing calculator!](https://www.desmos.com/calculator/2hm3b1fbdf)
```julia
p1 = 2
p2 = 2
k1 = KnotVector(-10:10)+p1*KnotVector(-10,10)
k2 = KnotVector(-10:10)+p2*KnotVector(-10,10)
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)

f(u1, u2) = [2u1 + sin(u1) + cos(u2) + u2 / 2, 3u2 + sin(u2) + sin(u1) / 2 + u1^2 / 6] / 5

a0 = fittingcontrolpoints(f, (P1, P2))
a = [a0[i1,i2][j] for i1 in 1:dim(P1), i2 in 1:dim(P2), j in 1:2]
M = BSplineManifold(a, (P1, P2))
save_png("fitting.png", M, unitlength=50, up=10, down=-10, left=-10, right=10)
```
![](docs/src/img/fitting_desmos.png)
![](docs/src/img/fitting.png)

If the knotvector span is too coarse, the approximation will be coarse.
```julia
p1 = 2
p2 = 2
k1 = KnotVector(-10:5:10)+p1*KnotVector(-10,10)
k2 = KnotVector(-10:5:10)+p2*KnotVector(-10,10)
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)

f(u1, u2) = [2u1 + sin(u1) + cos(u2) + u2 / 2, 3u2 + sin(u2) + sin(u1) / 2 + u1^2 / 6] / 5

a0 = fittingcontrolpoints(f, (P1, P2))
a = [a0[i1,i2][j] for i1 in 1:dim(P1), i2 in 1:dim(P2), j in 1:2]
M = BSplineManifold(a, (P1, P2))
save_png("fitting_coarse.png", M, unitlength=50, up=10, down=-10, left=-10, right=10)
```
![](docs/src/img/fitting_coarse.png)

### Draw smooth vector graphics
```julia
p = 3
k = KnotVector(range(-2π,2π,length=8))+p*KnotVector(-2π,2π)
P = BSplineSpace{p}(k)

f(u) = [u,sin(u)]

a0 = fittingcontrolpoints(f, (P,))
a = [a0[i1][j] for i1 in 1:dim(P), j in 1:2]
M = BSplineManifold(a, (P,))
save_svg("sine-curve.svg", M, unitlength=50, up=2, down=-2, left=-8, right=8)
save_svg("sine-curve_no-points.svg", M, unitlength=50, up=2, down=-2, left=-8, right=8, points=false)
```
![](docs/src/img/sine-curve.svg)
![](docs/src/img/sine-curve_no-points.svg)

This is useful when you edit graphs (or curves) with your favorite vector graphics editor.

![](docs/src/img/inkscape.png)

## References
If you use BasicBSpline.jl in your work, please consider citing it by

```bibtex
@misc{hyrodium:2020:BasicBSpline,
  title={BasicBSpline.jl: Basic operations for B-spline functions and related things with julia},
  author={Yuto Horikawa},
  year={2020},
  howpublished={\url{https://hyrodium.github.io/BasicBSpline.jl/stable/}},
  doi={10.5281/zenodo.4356869}
}
```
