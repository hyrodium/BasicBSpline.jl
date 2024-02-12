# Refinement

## Setup

```@example math_refinement
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots
```

## Example
### Define original manifold

```@example math_refinement
p = 2 # degree of polynomial
k = KnotVector(1:8) # knot vector
P = BSplineSpace{p}(k) # B-spline space
rand_a = [SVector(rand(), rand()) for i in 1:dim(P), j in 1:dim(P)]
a = [SVector(2*i-6.5, 2*j-6.5) for i in 1:dim(P), j in 1:dim(P)] + rand_a # random 
M = BSplineManifold(a,(P,P)) # Define B-spline manifold
gr()
plot(M)
savefig("refinement_2dim_original.png") # hide
nothing # hide
```

![](refinement_2dim_original.png)

### h-refinement
Insert additional knots to knot vector without changing the shape.

```@repl math_refinement
k₊ = (KnotVector([3.3,4.2]), KnotVector([3.8,3.2,5.3])) # additional knot vectors
M_h = refinement(M, k₊) # refinement of B-spline manifold
plot(M_h)
savefig("refinement_2dim_h.png") # hide
nothing # hide
```

![](refinement_2dim_h.png)

### p-refinement
Increase the polynomial degree of B-spline manifold without changing the shape.

```@repl math_refinement
p₊ = (Val(1), Val(2)) # additional degrees
M_p = refinement(M, p₊) # refinement of B-spline manifold
savefig("refinement_2dim_p.png") # hide
nothing # hide
```

![](refinement_2dim_p.png)
