# BasicBSpline

| **Documentation** | **Build Status** |
|:---:|:---:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hyrodium.github.io/BasicBSpline.jl/dev) | [![Build Status](https://travis-ci.com/hyrodium/BasicBSpline.jl.svg?branch=master)](https://travis-ci.com/hyrodium/BasicBSpline.jl) |

## Summary
This package provides basic (mathematical) operations for [B-spline](https://en.wikipedia.org/wiki/B-spline).

The package [Interpolations.jl](http://juliamath.github.io/Interpolations.jl/v0.12/) says:
> Currently this package's support is best for B-splines and also supports irregular grids. However, the API has been designed with intent to support more options. Pull-requests are more than welcome! It should be noted that the API may continue to evolve over time.

As mentioned before, this package treats mathematical aspect of B-spline, so the difference between these packages is a main purpose.
* If you are interested in Interpolations, [Interpolations.jl](http://juliamath.github.io/Interpolations.jl/v0.12/) would be helpful.
* If you would like to deal with *raw* B-spline functions, this package would be the best for you.  For example:
    * B-spline curve
    * B-spline surface
    * NURBS (Non-Uniform Rational B-Spline)
    * IGA (Isogeometric Analysis)

## Installation
```
(pkg) > add https://github.com/hyrodium/BasicBSpline.jl
```

## Example Images
###  Example of B-spline function


```julia
using BasicBSpline
using Plots
gr()

p = 3
k = Knots([0.00,1.50,2.50,5.50,8.00,9.00,9.50,10.0])
P0 = 𝒫(0,k)
P1 = 𝒫(1,k)
P2 = 𝒫(2,k)
P3 = 𝒫(3,k)
plot(
    plot([t->bsplinebasis(i,P0,t) for i in 1:dim(P0)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P1,t) for i in 1:dim(P1)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P2,t) for i in 1:dim(P2)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P3,t) for i in 1:dim(P3)], 0, 10, ylims=(0,1), legend=false),
    layout=(2,2),
)
```

![](docs/src/img/cover.png)

Try [interactive graph with desmos graphing calculator](https://www.desmos.com/calculator/ql6jqgdabs)!

![](docs/src/img/bsplinebasis.png)
