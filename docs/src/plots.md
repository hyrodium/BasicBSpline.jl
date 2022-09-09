# Plots.jl

`BasicBSpline.jl` has a dependency on [`RecipesBase.jl`](https://github.com/JuliaPlots/RecipesBase.jl).
This means, users can easily visalize instances defined in `BasicBSpline`.
In this section, we will provide some plottig examples.

```@setup plots
using BasicBSpline
using StaticArrays
using Plots; plotly()
```

## `BSplineSpace`

```@example plots
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P0 = BSplineSpace{0}(k) # 0th degree piecewise polynomial space
P1 = BSplineSpace{1}(k) # 1st degree piecewise polynomial space
P2 = BSplineSpace{2}(k) # 2nd degree piecewise polynomial space
P3 = BSplineSpace{3}(k) # 3rd degree piecewise polynomial space
plot(
    plot([t->bsplinebasis(P0,i,t) for i in 1:dim(P0)], 0, 10, ylims=(0,1), legend=false, title="0th polynomial degree"),
    plot([t->bsplinebasis(P1,i,t) for i in 1:dim(P1)], 0, 10, ylims=(0,1), legend=false, title="1st polynomial degree"),
    plot([t->bsplinebasis(P2,i,t) for i in 1:dim(P2)], 0, 10, ylims=(0,1), legend=false, title="2nd polynomial degree"),
    plot([t->bsplinebasis(P3,i,t) for i in 1:dim(P3)], 0, 10, ylims=(0,1), legend=false, title="3rd polynomial degree"),
)
savefig("plots-bsplinebasis-raw.html") # hide
```

```@raw html
<object type="text/html" data="../plots-bsplinebasis-raw.html" style="width:100%;height:420px;"></object>
```

```@example plots
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P0 = BSplineSpace{0}(k) # 0th degree piecewise polynomial space
P1 = BSplineSpace{1}(k) # 1st degree piecewise polynomial space
P2 = BSplineSpace{2}(k) # 2nd degree piecewise polynomial space
P3 = BSplineSpace{3}(k) # 3rd degree piecewise polynomial space
plot(
    plot(P0, ylims=(0,1), legend=false, title="0th polynomial degree"),
    plot(P1, ylims=(0,1), legend=false, title="1st polynomial degree"),
    plot(P2, ylims=(0,1), legend=false, title="2nd polynomial degree"),
    plot(P3, ylims=(0,1), legend=false, title="3rd polynomial degree"),
    layout=(2,2),
)
savefig("plots-bsplinebasis.html") # hide
```

```@raw html
<object type="text/html" data="../plots-bsplinebasis.html" style="width:100%;height:420px;"></object>
```

## `BSplineDerivativeSpace`

```@example plots
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P = BSplineSpace{3}(k)
plot(
    plot(BSplineDerivativeSpace{0}(P), label="0th derivative", color=:black),
    plot(BSplineDerivativeSpace{1}(P), label="1st derivative", color=:red),
    plot(BSplineDerivativeSpace{2}(P), label="2nd derivative", color=:green),
    plot(BSplineDerivativeSpace{3}(P), label="3rd derivative", color=:blue),
)
savefig("plots-bsplinebasisderivative.html") # hide
```

```@raw html
<object type="text/html" data="../plots-bsplinebasisderivative.html" style="width:100%;height:420px;"></object>
```

## `BSplineManifold`

### Cardioid (planar curve)
```@example plots
f(t) = SVector((1+cos(t))*cos(t),(1+cos(t))*sin(t))
p = 3
k = KnotVector(range(0,2π,15)) + p * KnotVector([0,2π]) + 2 * KnotVector([π])
P = BSplineSpace{p}(k)
a = fittingcontrolpoints(f,(P,))
M = BSplineManifold(a, (P,))

plot(M)
savefig("plots-cardioid.html") # hide
```

```@raw html
<object type="text/html" data="../plots-cardioid.html" style="width:100%;height:550px;"></object>
```

### Helix (spatial curve)
```@example plots
f(t) = SVector(cos(t),sin(t),t)
p = 3
k = KnotVector(range(0,6π,15)) + p * KnotVector([0,6π])
P = BSplineSpace{p}(k)
a = fittingcontrolpoints(f,(P,))
M = BSplineManifold(a, (P,))

plot(M)
savefig("plots-helix.html") # hide
```

```@raw html
<object type="text/html" data="../plots-helix.html" style="width:100%;height:550px;"></object>
```

### B-spline surface

```@example plots
p1 = 2
p2 = 3
k1 = KnotVector(1:10)
k2 = KnotVector(1:20)
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)
a = [SVector(i-j^2/20, j+i^2/10, sin((i+j)/2)+randn()) for i in 1:dim(P1), j in 1:dim(P2)]
M = BSplineManifold(a,(P1,P2))
plot(M)

savefig("plots-surface.html") # hide
```

```@raw html
<object type="text/html" data="../plots-surface.html" style="width:100%;height:550px;"></object>
```


## `RationalBSplineManifold`

```@example plots
k = KnotVector([0,0,0,1,1,1])
P = BSplineSpace{2}(k)
a = [SVector(1,0),SVector(1,1),SVector(0,1)]
w = [1,1/√2,1]
M = BSplineManifold(a,(P,))
R = RationalBSplineManifold(a,w,(P,))
ts = 0:0.01:2
plot(cospi.(ts),sinpi.(ts), label="circle")
plot!(M, label="B-spline curve")
plot!(R, label="Rational B-spline curve")

savefig("plots-arc.html") # hide
```

```@raw html
<object type="text/html" data="../plots-arc.html" style="width:100%;height:550px;"></object>
```
