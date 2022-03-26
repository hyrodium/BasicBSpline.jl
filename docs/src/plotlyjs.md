# PlotlyJS.jl

### sine curve

```@example
using BasicBSpline
using StaticArrays
using PlotlyJS
f(t) = SVector(t,sin(t))
p = 3
k = KnotVector(range(0,4π,15)) + p * KnotVector(0,4π)
P = BSplineSpace{p}(k)
contolpoints = fittingcontrolpoints(f,(P,))
M = BSplineManifold(contolpoints, (P,))

ts = range(0,4π,100)
xs_a = [a[1] for a in contolpoints]
ys_a = [a[2] for a in contolpoints]
xs_f = [M(t)[1] for t in ts]
ys_f = [M(t)[2] for t in ts]
fig = Plot(scatter(x=xs_a, y=ys_a, name="control points", mode="lines", line_color="blue"))
addtraces!(fig, scatter(x=xs_f, y=ys_f, name="B-spline curve", mode="lines", line_color="red"))
relayout!(fig, width=500, height=500)
savefig(fig,"a.html")
```

```@raw html
<iframe src="../a.html" width="500" height="500"></iframe>
```
