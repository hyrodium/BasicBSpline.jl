# B-spline manifold

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

## Multi-dimensional B-spline

!!! info "Thm.  Basis of tensor product of B-spline spaces"
    The tensor product of B-spline spaces ``\mathcal{P}[p^1,k^1]\otimes\mathcal{P}[p^2,k^2]`` linear space with the following basis.
    ```math
    \mathcal{P}[p^1,k^1]\otimes\mathcal{P}[p^2,k^2]
    = \operatorname*{span}_{i,j} (B_{(i,p^1,k^1)} \otimes B_{(j,p^2,k^2)})
    ```
    where the basis are defined as
    ```math
    (B_{(i,p^1,k^1)} \otimes B_{(j,p^2,k^2)})(t^1, t^2)
    = B_{(i,p^1,k^1)}(t^1) \cdot B_{(j,p^2,k^2)}(t^2)
    ```

Higher dimensional tensor products ``\mathcal{P}[p^1,k^1]\otimes\cdots\otimes\mathcal{P}[p^d,k^d]`` are defined similarly.

## B-spline manifold
B-spline manifold is a parametric representation of a shape.

!!! tip "Def.  B-spline manifold"
    For given ``d``-dimensional B-spline basis functions ``B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)}`` and given points ``\bm{a}_{i^1,\dots,i^d} \in \mathbb{R}^{\hat{d}}``, B-spline manifold is defined by following equality:
    ```math
    \bm{p}(t^1,\dots,t^d;\bm{a}_{i^1,\dots,i^d})
    =\sum_{i^1,\dots,i^d}(B_{(i^1,p^1,k^1)} \otimes \cdots \otimes B_{(i^d,p^d,k^d)})(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
    ```
    Where ``\bm{a}_{i^1,\dots,i^d}`` are called **control points**.

We will also write ``\bm{p}(t^1,\dots,t^d; \bm{a})``, ``\bm{p}(t^1,\dots,t^d)``, ``\bm{p}(t; \bm{a})`` or ``\bm{p}(t)`` for simplicity.

```@repl math
P1 = BSplineSpace{1}(KnotVector([0,0,1,1]))
P2 = BSplineSpace{1}(KnotVector([1,1,2,3,3]))
n1 = dim(P1) # 2
n2 = dim(P2) # 3
a = [SVector(i, j) for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂ array.
M = BSplineManifold(a, (P1, P2))
```

```@docs
BSplineManifold
```

```@docs
RationalBSplineManifold
```

### B-spline curve
```@example math
## 1-dim B-spline manifold
p = 2 # degree of polynomial
k = KnotVector(1:12) # knot vector
P = BSplineSpace{p}(k) # B-spline space
a = [SVector(i-5, 3*sin(i^2)) for i in 1:dim(P)] # control points
M = BSplineManifold(a, (P,)) # Define B-spline manifold
plot(M)
savefig("1dim-manifold.html") # hide
```

```@raw html
<object type="text/html" data="../1dim-manifold.html" style="width:100%;height:420px;"></object>
```

### B-spline surface
```@example math
## 2-dim B-spline manifold
p = 2 # degree of polynomial
k = KnotVector(1:8) # knot vector
P = BSplineSpace{p}(k) # B-spline space
rand_a = [SVector(rand(), rand(), rand()) for i in 1:dim(P), j in 1:dim(P)]
a = [SVector(2*i-6.5, 2*j-6.5, 0) for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
M = BSplineManifold(a,(P,P)) # Define B-spline manifold
plot(M)
savefig("2dim-manifold.html") # hide
```

```@raw html
<object type="text/html" data="../2dim-manifold.html" style="width:100%;height:420px;"></object>
```

![](2dim.png)

## Affine commutativity
!!! info "Thm.  Affine commutativity"
    If ``T`` is a affine transform ``\mathbb{R}^d\to\mathbb{R}^d``, then the following equality holds.
    ```math
    T(\bm{p}(t; \bm{a}))
    =\bm{p}(t; T(\bm{a}))
    ```

