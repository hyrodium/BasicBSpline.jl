# PlotlyJS.jl

### Cardioid

```@example
using BasicBSpline
using StaticArrays
using PlotlyJS
f(t) = SVector((1+cos(t))*cos(t),(1+cos(t))*sin(t))
p = 3
k = KnotVector(range(0,2π,15)) + p * KnotVector(0,2π) + 2 * KnotVector(π)
P = BSplineSpace{p}(k)
contolpoints = fittingcontrolpoints(f,(P,))
M = BSplineManifold(contolpoints, (P,))

ts = range(0,2π,250)
xs_a = [a[1] for a in contolpoints]
ys_a = [a[2] for a in contolpoints]
xs_f = [M(t)[1] for t in ts]
ys_f = [M(t)[2] for t in ts]
fig = Plot(scatter(x=xs_a, y=ys_a, name="control points", line_color="blue", marker_size=8))
addtraces!(fig, scatter(x=xs_f, y=ys_f, name="B-spline curve", mode="lines", line_color="red"))
relayout!(fig, width=500, height=500)
savefig(fig,"cardioid.html")
nothing # hide
```

```@raw html
<object type="text/html" data="../cardioid.html" style="width:100%;height:550px;"></object>
```
