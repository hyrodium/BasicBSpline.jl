# Interpolations

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

# Example inputs
xs = [1, 2, 3, 4, 6, 7]
fs = [1.3, 1.5, 2, 2.1, 1.9, 1.3]
f = interpolate(xs,fs)

# Plot
scatter(xs, fs)
plot!(t->f(t))
savefig("interpolation_cubic.html") # hide
```

```@raw html
<object type="text/html" data="../interpolation_cubic.html" style="width:100%;height:420px;"></object>
```

## Interpolation with linear B-spline

```@example interpolation
function interpolate_linear(xs::AbstractVector, fs::AbstractVector{T}) where T
    # Linear open B-spline space
    p = 1
    k = KnotVector(xs) + KnotVector(xs[1],xs[end])
    P = BSplineSpace{p}(k)

    # dimensions
    m = length(xs)
    n = dim(P)

    # Compute the interpolant function (1-dim B-spline manifold)
    return BSplineManifold(fs, (P,))
end

# Example inputs
xs = [1, 2, 3, 4, 6, 7]
fs = [1.3, 1.5, 2, 2.1, 1.9, 1.3]

f = interpolate_linear(xs,fs)

# Plot
scatter(xs, fs)
plot!(t->f(t))
savefig("interpolation_linear.html") # hide
```

```@raw html
<object type="text/html" data="../interpolation_linear.html" style="width:100%;height:420px;"></object>
```

## Interpolation with periodic B-spline

```@example interpolation
function interpolate_periodic(xs::AbstractVector, fs::AbstractVector, ::Val{p}) where p
    # Closed B-spline space, any polynomial degrees can be accepted
    n = length(xs) - 1
    period = xs[end]-xs[begin]
    k = KnotVector(vcat(
        xs[end-p:end-1] .- period,
        xs,
        xs[begin+1:begin+p] .+ period
    ))
    P = BSplineSpace{p}(k)
    A = [bsplinebasis(P,j,xs[i]) for i in 1:n, j in 1:n]
    for i in 1:p-1, j in 1:i
        A[n+i-p+1,j] += bsplinebasis(P,j+n,xs[i+n-p+1])
    end
    b = A \ fs[begin:end-1]
    # Compute the interpolant function (1-dim B-spline manifold)
    return BSplineManifold(vcat(b,b[1:p]), (P,))
end

# Example inputs
xs = [1, 2, 3, 4, 6, 7]
fs = [1.3, 1.5, 2, 2.1, 1.9, 1.3] # fs[1] == fs[end]

f = interpolate_periodic(xs,fs,Val(2))

# Plot
scatter(xs, fs)
plot!(t->f(mod(t-1,6)+1),1,14)
plot!(t->f(t))
savefig("interpolation_periodic.html") # hide
```

```@raw html
<object type="text/html" data="../interpolation_periodic.html" style="width:100%;height:420px;"></object>
```

Note that the periodic interpolation supports any degree of polynomial.

```@example interpolation
xs = 2π*rand(5)
sort!(push!(xs, 0, 2π))
fs = sin.(xs)
scatter(xs, fs)
f1 = interpolate_periodic(xs,fs,Val(1))
f2 = interpolate_periodic(xs,fs,Val(2))
f3 = interpolate_periodic(xs,fs,Val(3))
f4 = interpolate_periodic(xs,fs,Val(4))
f5 = interpolate_periodic(xs,fs,Val(5))
plot!(t->f1(t), label="polynomial degree 1")
plot!(t->f2(t), label="polynomial degree 2")
plot!(t->f3(t), label="polynomial degree 3")
plot!(t->f4(t), label="polynomial degree 4")
plot!(t->f5(t), label="polynomial degree 5")
savefig("interpolation_periodic_sin.html") # hide
```

```@raw html
<object type="text/html" data="../interpolation_periodic_sin.html" style="width:100%;height:420px;"></object>
```
