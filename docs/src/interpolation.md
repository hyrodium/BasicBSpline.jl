# Interpolation

Currently, BasicBSpline.jl doesn't have APIs for interpolations, but it is not hard to implement the basic interpolation algorithms.

```@setup interpolation
using BasicBSpline
using IntervalSets
using Plots; plotly()
```

## Interpolation with cubic B-spline

```@example interpolation
function interpolate(xs::AbstractVector, fs::AbstractVector{T}) where T
    # Cubic open B-spline space
    p = 3
    k = KnotVector(xs) + KnotVector(xs[1],xs[end]) * p
    P = BSplineSpace{p}(k)

    # dimensions
    m = length(xs)
    n = dim(P)

    # The interpolant function has a f''=0 property at bounds.
    ddP = BSplineDerivativeSpace{2}(P)
    dda = [bsplinebasis(ddP,j,xs[1]) for j in 1:n]
    ddb = [bsplinebasis(ddP,j,xs[m]) for j in 1:n]

    # Compute the interpolant function (1-dim B-spline manifold)
    M = [bsplinebasis(P,j,xs[i]) for i in 1:m, j in 1:n]
    M = vcat(dda', M, ddb')
    y = vcat(zero(T), fs, zero(T))
    return BSplineManifold(M\y, (P,))
end

## Example inputs
xs = [1, 2, 3, 4, 6, 7]
fs = [1.3, 1.5, 2, 2.1, 1.9, 1.3]
f = interpolate(xs,fs)

## Plot
scatter(xs, fs)
plot!(t->f(t))
savefig("interpolation_cubic.html") # hide
```

```@raw html
<object type="text/html" data="../interpolation_cubic.html" style="width:100%;height:420px;"></object>
```
