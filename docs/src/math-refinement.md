# Refinement

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots
```

```@docs
refinement
```

```@repl math
p = 2 # degree of polynomial
k = KnotVector(1:8) # knot vector
P = BSplineSpace{p}(k) # B-spline space
rand_a = [SVector(rand(), rand()) for i in 1:dim(P), j in 1:dim(P)]
a = [SVector(2*i-6.5, 2*j-6.5) for i in 1:dim(P), j in 1:dim(P)] + rand_a # random 
M = BSplineManifold(a,(P,P)) # Define B-spline manifold
```

## h-refinemnet
Insert additional knots to knot vector.

```@repl math
k₊=(KnotVector(3.3,4.2),KnotVector(3.8,3.2,5.3)) # additional knotvectors
M_h = refinement(M,k₊=k₊) # refinement of B-spline manifold
save_png("2dim_h-refinement.png", M_h) # save image
```
![](2dim_h-refinement.png)

Note that this shape and the last shape are identical.


## p-refinemnet
Increase the polynomial degree of B-spline manifold.

```@repl math
p₊=(1,2) # additional degrees
M_p = refinement(M,p₊=p₊) # refinement of B-spline manifold
save_png("2dim_p-refinement.png", M_p) # save image
```
![](2dim_p-refinement.png)

Note that this shape and the last shape are identical.
