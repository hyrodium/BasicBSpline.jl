# Geometric modeling

## Load packages

```@example geometricmodeling
using BasicBSpline
using StaticArrays
using Plots
using LinearAlgebra
plotly()
```

## Arc
```@example geometricmodeling
p = 2
k = KnotVector([0,0,0,1,1,1])
P = BSplineSpace{p}(k)
t = 1  # angle in radians
a = [SVector(1,0), SVector(1,tan(t/2)), SVector(cos(t),sin(t))]
w = [1,cos(t/2),1]
M = RationalBSplineManifold(a,w,P)
plot(M, xlims=(0,1.1), ylims=(0,1.1), aspectratio=1)
savefig("geometricmodeling-arc.html") # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-arc.html" style="width:100%;height:420px;"></object>
```

## Circle
```@example geometricmodeling
p = 2
k = KnotVector([0,0,0,1,1,2,2,3,3,4,4,4])
P = BSplineSpace{p}(k)
a = [
    SVector(1,0),
    SVector(1,1),
    SVector(0,1),
    SVector(-1,1),
    SVector(-1,0),
    SVector(-1,-1),
    SVector(0,-1),
    SVector(1,-1),
    SVector(1,0)
]
w = [1,1/√2,1,1/√2,1,1/√2,1,1/√2,1]
M = RationalBSplineManifold(a,w,P)
plot(M, xlims=(-1.2,1.2), ylims=(-1.2,1.2), aspectratio=1)
savefig("geometricmodeling-circle.html") # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-circle.html" style="width:100%;height:420px;"></object>
```

## Circle
```@example geometricmodeling
R = 3
r = 1

a1 = (R+r)*a
a5 = (R-r)*a
a2 = [p+r*SVector(0,0,1) for p in a1]
a3 = [p+r*SVector(0,0,1) for p in R*a]
a4 = [p+r*SVector(0,0,1) for p in a5]
a6 = [p-r*SVector(0,0,1) for p in a5]
a7 = [p-r*SVector(0,0,1) for p in R*a]
a8 = [p-r*SVector(0,0,1) for p in a1]
a9 = a1

A = hcat(a1,a2,a3,a4,a5,a6,a7,a8,a9)
M = RationalBSplineManifold(A,w*w',P,P)
plot(M)
savefig("geometricmodeling-torus.html") # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-torus.html" style="width:100%;height:420px;"></object>
```

## Paraboloid
```@example geometricmodeling
p = 2
k = KnotVector([-1,-1,-1,1,1,1])
P = BSplineSpace{p}(k)
a = [SVector(i,j,2i^2+2j^2-2) for i in -1:1, j in -1:1]
M = BSplineManifold(a,P,P)
plot(M)
savefig("geometricmodeling-paraboloid.html") # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-paraboloid.html" style="width:100%;height:420px;"></object>
```

## Hyperbolic paraboloid
```@example geometricmodeling
a = [SVector(i,j,2i^2-2j^2) for i in -1:1, j in -1:1]
M = BSplineManifold(a,P,P)
plot(M)
savefig("geometricmodeling-hyperbolicparaboloid.html") # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-hyperbolicparaboloid.html" style="width:100%;height:420px;"></object>
```
