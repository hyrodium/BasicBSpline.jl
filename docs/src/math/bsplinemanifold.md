# B-spline manifold

## Setup
```@example math_bsplinemanifold
using BasicBSpline
using StaticArrays
using StaticArrays
using Plots; plotly()
```

## Multi-dimensional B-spline

!!! info "Thm.  Basis of tensor product of B-spline spaces"
    The tensor product of B-spline spaces ``\mathcal{P}[p^1,k^1]\otimes\mathcal{P}[p^2,k^2]`` is a linear space with the following basis.
    ```math
    \mathcal{P}[p^1,k^1]\otimes\mathcal{P}[p^2,k^2]
    = \operatorname*{span}_{i,j} (B_{(i,p^1,k^1)} \otimes B_{(j,p^2,k^2)})
    ```
    where the basis are defined as
    ```math
    (B_{(i,p^1,k^1)} \otimes B_{(j,p^2,k^2)})(t^1, t^2)
    = B_{(i,p^1,k^1)}(t^1) \cdot B_{(j,p^2,k^2)}(t^2)
    ```

The next plot shows ``B_{(3,p^1,k^1)} \otimes B_{(4,p^2,k^2)}`` basis function.

```@example math_bsplinemanifold
# Define shape
k1 = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
k2 = knotvector"31 2  121"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{2}(k2)

# Choose indices
i1 = 3
i2 = 4

# Visualize basis functions
plot(P1; plane=:xz, label="P1", color=:red)
plot!(P2; plane=:yz, label="P2", color=:green)
plot!(k1; plane=:xz, label="k1", color=:red, markersize=2)
plot!(k2; plane=:yz, label="k2", color=:green, markersize=2)
xs = range(0,10,length=100)
ys = range(1,9,length=100)
surface!(xs, ys, bsplinebasis.(P1,i1,xs') .* bsplinebasis.(P2,i2,ys))
savefig("2dim-tensor-product-bspline.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../2dim-tensor-product-bspline.html" style="width:100%;height:420px;"></object>
```

Higher dimensional tensor products ``\mathcal{P}[p^1,k^1]\otimes\cdots\otimes\mathcal{P}[p^d,k^d]`` are defined similarly.

## B-spline manifold
B-spline manifold is a parametric representation of a shape.

!!! tip "Def.  B-spline manifold"
    For given ``d``-dimensional B-spline basis functions ``B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)}`` and given points ``\bm{a}_{i^1 \dots i^d} \in V``, B-spline manifold is defined by the following parametrization:
    ```math
    \bm{p}(t^1,\dots,t^d;\bm{a}_{i^1 \dots i^d})
    =\sum_{i^1,\dots,i^d}(B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)})(t^1,\dots,t^d) \bm{a}_{i^1 \dots i^d}
    ```
    Where ``\bm{a}_{i^1 \dots i^d}`` are called **control points**.

We will also write ``\bm{p}(t^1,\dots,t^d; \bm{a})``, ``\bm{p}(t^1,\dots,t^d)``, ``\bm{p}(t; \bm{a})`` or ``\bm{p}(t)`` for simplicity.

Note that the [`BSplineManifold`](@ref) objects are callable, and the arguments will be checked if it fits in the domain of [`BSplineSpace`](@ref).

```@example math_bsplinemanifold
P = BSplineSpace{2}(KnotVector([0,0,0,1,1,1]))
a = [SVector(1,0), SVector(1,1), SVector(0,1)]  # `length(a) == dim(P)`
M = BSplineManifold(a, P)
M(0.4)  # Calculate `sum(a[i]*bsplinebasis(P, i, 0.4) for i in 1:dim(P))`
```

If you need extension of [`BSplineManifold`](@ref) or don't need the arguments check, you can call [`unbounded_mapping`](@ref).

```@repl math_bsplinemanifold
M(0.4)
unbounded_mapping(M, 0.4)
M(1.2)
unbounded_mapping(M, 1.2)
```

`unbounded_mapping(M,t...)` is a little bit faster than `M(t...)` because it does not check the domain.

### B-spline curve
```@example math_bsplinemanifold
## 1-dim B-spline manifold
p = 2 # degree of polynomial
k = KnotVector(1:12) # knot vector
P = BSplineSpace{p}(k) # B-spline space
a = [SVector(i-5, 3*sin(i^2)) for i in 1:dim(P)] # control points
M = BSplineManifold(a, P) # Define B-spline manifold
plot(M)
savefig("1dim-manifold.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../1dim-manifold.html" style="width:100%;height:420px;"></object>
```

### B-spline surface
```@example math_bsplinemanifold
## 2-dim B-spline manifold
p = 2 # degree of polynomial
k = KnotVector(1:8) # knot vector
P = BSplineSpace{p}(k) # B-spline space
rand_a = [SVector(rand(), rand(), rand()) for i in 1:dim(P), j in 1:dim(P)]
a = [SVector(2*i-6.5, 2*j-6.5, 0) for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
M = BSplineManifold(a,(P,P)) # Define B-spline manifold
plot(M)
savefig("2dim-manifold.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../2dim-manifold.html" style="width:100%;height:420px;"></object>
```

## Affine commutativity
!!! info "Thm.  Affine commutativity"
    Let ``T`` be a affine transform ``V \to W``, then the following equality holds.
    ```math
    T(\bm{p}(t; \bm{a}))
    =\bm{p}(t; T(\bm{a}))
    ```

## Fixing arguments (currying)
Just like fixing first index such as `A[3,:]::Array{1}` for matrix `A::Array{2}`, fixing first argument `M(4.3,:)` will create `BSplineManifold{1}` for B-spline surface `M::BSplineManifold{2}`.

```@repl math_bsplinemanifold
p = 2;
k = KnotVector(1:8);
P = BSplineSpace{p}(k);
a = [SVector(5i+5j+rand(), 5i-5j+rand(), rand()) for i in 1:dim(P), j in 1:dim(P)];
M = BSplineManifold(a,(P,P));
M isa BSplineManifold{2}
M(:,:) isa BSplineManifold{2}
M(4.3,:) isa BSplineManifold{1}  # Fix first argument
```

```@example math_bsplinemanifold
plot(M)
plot!(M(4.3,:), linewidth = 5, color=:cyan)
plot!(M(4.4,:), linewidth = 5, color=:red)
plot!(M(:,5.2), linewidth = 5, color=:green)
savefig("2dim-manifold-currying.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../2dim-manifold-currying.html" style="width:100%;height:420px;"></object>
```
