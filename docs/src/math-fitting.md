# Fitting with B-spline manifold

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

Fitting with least squares method.

```@docs
fittingcontrolpoints
```

[Try on Desmos graphing graphing calculator!](https://www.desmos.com/calculator/2hm3b1fbdf)
```@repl math
p1 = 2
p2 = 2
k1 = KnotVector(-10:10)+p1*KnotVector(-10,10)
k2 = KnotVector(-10:10)+p2*KnotVector(-10,10)
P1 = BSplineSpace{p1}(k1)
P2 = BSplineSpace{p2}(k2)

f(u1, u2) = SVector(2u1 + sin(u1) + cos(u2) + u2 / 2, 3u2 + sin(u2) + sin(u1) / 2 + u1^2 / 6) / 5

a = fittingcontrolpoints(f, (P1, P2))
M = BSplineManifold(a, (P1, P2))
save_png("fitting.png", M, unitlength=50, xlims=(-10,10), ylims=(-10,10))
```
![](img/fitting_desmos.png)
![](fitting.png)
