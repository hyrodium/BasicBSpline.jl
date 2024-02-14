# BasicBSpline.jl

Basic (mathematical) operations for B-spline functions and related things with Julia.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/dev)
[![Build Status](https://github.com/hyrodium/BasicBSpline.jl/workflows/CI/badge.svg)](https://github.com/hyrodium/BasicBSpline.jl/actions)
[![Coverage](https://codecov.io/gh/hyrodium/BasicBSpline.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hyrodium/BasicBSpline.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/258791290.svg)](https://zenodo.org/badge/latestdoi/258791290)
[![BasicBSpline Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/BasicBSpline)](https://pkgs.genieframework.com?packages=BasicBSpline)

![](docs/src/img/BasicBSplineLogo.png)

## Summary
This package provides basic mathematical operations for [B-spline](https://en.wikipedia.org/wiki/B-spline).

* B-spline basis function
* Some operations for knot vector
* Some operations for B-spline space (piecewise polynomial space)
* B-spline manifold (includes curve, surface and solid)
* Refinement algorithm for B-spline manifold
* Fitting control points for a given function

## Comparison to other B-spline packages
There are several Julia packages for B-spline, and this package distinguishes itself with the following key benefits:

* Supports all degrees of polynomials.
* Includes a refinement algorithm for B-spline manifolds.
* Delivers high-speed performance.
* Is mathematically oriented.
* Provides a fitting algorithm using least squares. (via [BasicBSplineFitting.jl](https://github.com/hyrodium/BasicBSplineFitting.jl))
* Offers exact SVG export feature. (via [BasicBSplineExporter.jl](https://github.com/hyrodium/BasicBSplineExporter.jl))

If you have any thoughts, please comment in:

* [Discourse post about BasicBSpline.jl](https://discourse.julialang.org/t//96323)
* [JuliaPackageComparisons](https://juliapackagecomparisons.github.io/pages/bspline/)

## Installation
Install this package via Julia REPL's package mode.

```
]add BasicBSpline
```

## Quick start
### B-spline basis function

The value of B-spline basis function $B_{(i,p,k)}$ can be obtained with `bsplinebasis₊₀`.

$$
\begin{aligned}
{B}_{(i,p,k)}(t)
&=
\frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
+\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}\le t< k_{i+1})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
$$

```julia
julia> using BasicBSpline

julia> P3 = BSplineSpace{3}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))
BSplineSpace{3, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> bsplinebasis₊₀(P3, 2, 7.5)
0.13786213786213783
```

BasicBSpline.jl has many recipes based on [RecipesBase.jl](https://docs.juliaplots.org/dev/RecipesBase/), and `BSplineSpace` object can be visualized with its basis functions.
([Try B-spline basis functions in Desmos](https://www.desmos.com/calculator/ql6jqgdabs))

```julia
using BasicBSpline
using Plots

k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P0 = BSplineSpace{0}(k) # 0th degree piecewise polynomial space
P1 = BSplineSpace{1}(k) # 1st degree piecewise polynomial space
P2 = BSplineSpace{2}(k) # 2nd degree piecewise polynomial space
P3 = BSplineSpace{3}(k) # 3rd degree piecewise polynomial space

gr()
plot(
    plot(P0; ylims=(0,1), label="P0"),
    plot(P1; ylims=(0,1), label="P1"),
    plot(P2; ylims=(0,1), label="P2"),
    plot(P3; ylims=(0,1), label="P3"),
    layout=(2,2),
)
```
![](docs/src/img/cover.png)

You can visualize the differentiability of B-spline basis function. See [Differentiability and knot duplications](https://hyrodium.github.io/BasicBSpline.jl/dev/math/bsplinebasis/#Differentiability-and-knot-duplications) for details.

https://github.com/hyrodium/BasicBSpline.jl/assets/7488140/034cf6d0-62ea-44e0-a00f-d307e6aad0fe

### B-spline manifold
```julia
using BasicBSpline
using StaticArrays
using Plots

## 1-dim B-spline manifold
p = 2 # degree of polynomial
k = KnotVector(1:12) # knot vector
P = BSplineSpace{p}(k) # B-spline space
a = [SVector(i-5, 3*sin(i^2)) for i in 1:dim(P)] # control points
M = BSplineManifold(a, P) # Define B-spline manifold
gr(); plot(M)
```
![](docs/src/img/bspline_curve.png)

### Rational B-spline manifold (NURBS)

```julia
using BasicBSpline
using LinearAlgebra
using StaticArrays
using Plots
plotly()

R1 = 3  # major radius of torus
R2 = 1  # minor radius of torus

p = 2
k = KnotVector([0,0,0,1,1,2,2,3,3,4,4,4])
P = BSplineSpace{p}(k)
a = [normalize(SVector(cosd(t), sind(t)), Inf) for t in 0:45:360]
w = [ifelse(isodd(i), √2, 1) for i in 1:9]

a0 = push.(a, 0)
a1 = (R1+R2)*a0
a5 = (R1-R2)*a0
a2 = [p+R2*SVector(0,0,1) for p in a1]
a3 = [p+R2*SVector(0,0,1) for p in R1*a0]
a4 = [p+R2*SVector(0,0,1) for p in a5]
a6 = [p-R2*SVector(0,0,1) for p in a5]
a7 = [p-R2*SVector(0,0,1) for p in R1*a0]
a8 = [p-R2*SVector(0,0,1) for p in a1]

M = RationalBSplineManifold(hcat(a1,a2,a3,a4,a5,a6,a7,a8,a1), w*w', P, P)
plot(M; controlpoints=(markersize=2,))
```
![](docs/src/img/rational_bspline_surface_plotly.png)

### Refinement
#### h-refinement (knot insertion)
Insert additional knots to knot vectors without changing the shape.

```julia
k₊ = (KnotVector([3.1, 3.2, 3.3]), KnotVector([0.5, 0.8, 0.9])) # additional knot vectors
M_h = refinement(M, k₊) # refinement of B-spline manifold
plot(M_h; controlpoints=(markersize=2,))
```
![](docs/src/img/rational_bspline_surface_href_plotly.png)

#### p-refinement (degree elevation)
Increase the polynomial degrees of B-spline manifold without changing the shape.

```julia
p₊ = (Val(1), Val(2)) # additional degrees
M_p = refinement(M, p₊) # refinement of B-spline manifold
plot(M_p; controlpoints=(markersize=2,))
```
![](docs/src/img/rational_bspline_surface_pref_plotly.png)

### Fitting B-spline manifold
The next example shows the fitting for [the following graph on Desmos graphing calculator](https://www.desmos.com/calculator/2hm3b1fbdf)!

![](docs/src/img/fitting_desmos.png)

```julia
using BasicBSplineFitting

p1 = 2
p2 = 2
k1 = KnotVector(-10:10)+p1*KnotVector([-10,10])
k2 = KnotVector(-10:10)+p2*KnotVector([-10,10])
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)

f(u1, u2) = SVector(2u1 + sin(u1) + cos(u2) + u2 / 2, 3u2 + sin(u2) + sin(u1) / 2 + u1^2 / 6) / 5

a = fittingcontrolpoints(f, (P1, P2))
M = BSplineManifold(a, (P1, P2))
gr()
plot(M; aspectratio=1)
```
![](docs/src/img/fitting.png)

If the knot vector span is too coarse, the approximation will be coarse.
```julia
p1 = 2
p2 = 2
k1 = KnotVector(-10:5:10)+p1*KnotVector([-10,10])
k2 = KnotVector(-10:5:10)+p2*KnotVector([-10,10])
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)

f(u1, u2) = SVector(2u1 + sin(u1) + cos(u2) + u2 / 2, 3u2 + sin(u2) + sin(u1) / 2 + u1^2 / 6) / 5

a = fittingcontrolpoints(f, (P1, P2))
M = BSplineManifold(a, (P1, P2))
plot(M; aspectratio=1)
```
![](docs/src/img/fitting_coarse.png)

### Draw smooth vector graphics
```julia
using BasicBSpline
using BasicBSplineFitting
using BasicBSplineExporter
p = 3
k = KnotVector(range(-2π,2π,length=8))+p*KnotVector([-2π,2π])
P = BSplineSpace{p}(k)

f(u) = SVector(u,sin(u))

a = fittingcontrolpoints(f, P)
M = BSplineManifold(a, P)
save_svg("sine-curve.svg", M, unitlength=50, xlims=(-8,8), ylims=(-2,2))
save_svg("sine-curve_no-points.svg", M, unitlength=50, xlims=(-8,8), ylims=(-2,2), points=false)
```
![](docs/src/img/sine-curve.svg)
![](docs/src/img/sine-curve_no-points.svg)

This is useful when you edit graphs (or curves) with your favorite vector graphics editor.

![](docs/src/img/inkscape.png)

See [Plotting smooth graphs with Julia](https://forem.julialang.org/hyrodium/plotting-smooth-graphs-with-julia-6mj) for more tutorials.
