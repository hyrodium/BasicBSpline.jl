# BasicBSpline.jl

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
Install this package

```julia
(pkg)> add https://github.com/hyrodium/BasicBSpline.jl
```

To export graphics, use [ExportNURBS.jl](https://github.com/hyrodium/ExportNURBS.jl).

```julia
(pkg)> add https://github.com/hyrodium/ExportNURBS.jl
```

## Example
### B-spline function

```julia
using BasicBSpline
using Plots
gr()

k = Knots([0.00,1.50,2.50,5.50,8.00,9.00,9.50,10.0])
P0 = BSplineSpace(0,k) # 0th degree piecewise polynomial space
P1 = BSplineSpace(1,k) # 1st degree piecewise polynomial space
P2 = BSplineSpace(2,k) # 2nd degree piecewise polynomial space
P3 = BSplineSpace(3,k) # 3rd degree piecewise polynomial space
plot(
    plot([t->bsplinebasis(i,P0,t) for i in 1:dim(P0)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P1,t) for i in 1:dim(P1)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P2,t) for i in 1:dim(P2)], 0, 10, ylims=(0,1), legend=false),
    plot([t->bsplinebasis(i,P3,t) for i in 1:dim(P3)], 0, 10, ylims=(0,1), legend=false),
    layout=(2,2),
)
```

![](img/cover.png)

Try [interactive graph with Desmos graphing calculator](https://www.desmos.com/calculator/ql6jqgdabs)!

### B-spline manifold
```julia
using BasicBSpline
using ExportNURBS

p = 2 # degree of polynomial
k = Knots(1:8) # knot vector
P = BSplineSpace(p,k) # B-spline space
rand_a = [rand(2) for i in 1:dim(P), j in 1:dim(P)]
a = [[2*i-6.5,2*j-6.5] for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
M = BSplineManifold([P,P],a) # Define B-spline manifold
save_png("docs/src/img/2dim.png", M) # save image
```
![](img/2dim.png)

### Refinement
```julia
k₊=[Knots(3.3,4.2),Knots(3.8,3.2,5.3)] # additional knots
M′ = refinement(M,k₊=k₊) # refinement of B-spline manifold
save_png("docs/src/img/2dim_refinement.png", M′) # save image
```
![](img/2dim_refinement.png)

Note that this shape and the last shape are identical.

### Fitting B-spline manifold
[Try on Desmos graphing graphing calculator!](https://www.desmos.com/calculator/2hm3b1fbdf)
```julia
p1 = 2
p2 = 2
k1 = Knots(-10:10)+p1*Knots(-10,10)
k2 = Knots(-10:10)+p2*Knots(-10,10)
P1 = FastBSplineSpace(p1, k1)
P2 = FastBSplineSpace(p2, k2)

f(u) = [2u[1]+sin(u[1])+cos(u[2])+u[2]/2, 3u[2]+sin(u[2])+sin(u[1])/2+u[1]^2/6]/5

a = fittingcontrolpoints(f, [P1,P2])
M = BSplineManifold([P1,P2],a)
save_png("docs/src/img/fitting.png", M, unitlength=50, up=10, down=-10, left=-10, right=10)
```
![](img/fitting_desmos.png)
![](img/fitting.png)

If the knots span is too coarse, the approximation will be coarse.
```julia
p1 = 2
p2 = 2
k1 = Knots(-10:5:10)+p1*Knots(-10,10)
k2 = Knots(-10:5:10)+p2*Knots(-10,10)
P1 = FastBSplineSpace(p1, k1)
P2 = FastBSplineSpace(p2, k2)

f(u) = [2u[1]+sin(u[1])+cos(u[2])+u[2]/2, 3u[2]+sin(u[2])+sin(u[1])/2+u[1]^2/6]/5

a = fittingcontrolpoints(f, [P1,P2])
M = BSplineManifold([P1,P2],a)
save_png("docs/src/img/fitting_coarse.png", M, unitlength=50, up=10, down=-10, left=-10, right=10)
```
![](img/fitting_coarse.png)

### Draw smooth vector graphics
```julia
p = 3
k = Knots(range(-2π,2π,length=8))+p*Knots(-2π,2π)
P = FastBSplineSpace(p, k)

f(u) = [u[1],sin(u[1])]

a = fittingcontrolpoints(f, [P])
M = BSplineManifold([P],a)
save_svg("docs/src/img/sine_curve.svg", M, unitlength=50, up=2, down=-2, left=-8, right=8)
```
![](img/sine_curve.svg)

This feature is useful when you use vector graphics editor.
