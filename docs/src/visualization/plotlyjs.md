# PlotlyJS.jl

## Cardioid

```@example
using BasicBSpline
using BasicBSplineFitting
using StaticArrays
using PlotlyJS
f(t) = SVector((1+cos(t))*cos(t),(1+cos(t))*sin(t))
p = 3
k = KnotVector(range(0,2π,15)) + p * KnotVector([0,2π]) + 2 * KnotVector([π])
P = BSplineSpace{p}(k)
a = fittingcontrolpoints(f, P)
M = BSplineManifold(a,  P)

ts = range(0,2π,250)
xs_a = getindex.(a,1)
ys_a = getindex.(a,2)
xs_f = getindex.(M.(ts),1)
ys_f = getindex.(M.(ts),2)
fig = Plot(scatter(x=xs_a, y=ys_a, name="control points", line_color="blue", marker_size=8))
addtraces!(fig, scatter(x=xs_f, y=ys_f, name="B-spline curve", mode="lines", line_color="red"))
relayout!(fig, width=500, height=500)
savefig(fig,"cardioid.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../cardioid.html" style="width:100%;height:550px;"></object>
```

## Helix
```@example
using BasicBSpline
using BasicBSplineFitting
using StaticArrays
using PlotlyJS
f(t) = SVector(cos(t),sin(t),t)
p = 3
k = KnotVector(range(0,6π,15)) + p * KnotVector([0,6π])
P = BSplineSpace{p}(k)
a = fittingcontrolpoints(f, P)
M = BSplineManifold(a,  P)

ts = range(0,6π,250)
xs_a = getindex.(a,1)
ys_a = getindex.(a,2)
zs_a = getindex.(a,3)
xs_f = getindex.(M.(ts),1)
ys_f = getindex.(M.(ts),2)
zs_f = getindex.(M.(ts),3)
fig = Plot(scatter3d(x=xs_a, y=ys_a, z=zs_a, name="control points", line_color="blue", marker_size=8))
addtraces!(fig, scatter3d(x=xs_f, y=ys_f, z=zs_f, name="B-spline curve", mode="lines", line_color="red"))
relayout!(fig, width=500, height=500)
savefig(fig,"helix.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../helix.html" style="width:100%;height:550px;"></object>
```
