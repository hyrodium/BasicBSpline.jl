# Geometric modeling

## Setup

```@example geometricmodeling
using BasicBSpline
using StaticArrays
using Plots
using LinearAlgebra
gr()
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
savefig("geometricmodeling-arc.png") # hide
nothing # hide
```

![](geometricmodeling-arc.png)

## Circle
```@example geometricmodeling
p = 2
k = KnotVector([0,0,0,1,1,2,2,3,3,4,4,4])
P = BSplineSpace{p}(k)
a = [normalize(SVector(cosd(t), sind(t)), Inf) for t in 0:45:360]
w = [ifelse(isodd(i), √2, 1) for i in 1:9]
M = RationalBSplineManifold(a,w,P)
plot(M, xlims=(-1.2,1.2), ylims=(-1.2,1.2), aspectratio=1)
savefig("geometricmodeling-circle.png") # hide
nothing # hide
```

![](geometricmodeling-circle.png)

## Torus
```@example geometricmodeling
plotly()
R1 = 3
R2 = 1

A = push.(a, 0)

a1 = (R1+R2)*A
a5 = (R1-R2)*A
a2 = [p+R2*SVector(0,0,1) for p in a1]
a3 = [p+R2*SVector(0,0,1) for p in R1*A]
a4 = [p+R2*SVector(0,0,1) for p in a5]
a6 = [p-R2*SVector(0,0,1) for p in a5]
a7 = [p-R2*SVector(0,0,1) for p in R1*A]
a8 = [p-R2*SVector(0,0,1) for p in a1]

a = hcat(a1,a2,a3,a4,a5,a6,a7,a8,a1)
M = RationalBSplineManifold(a,w*w',P,P)
plot(M)
savefig("geometricmodeling-torus.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-torus.html" style="width:100%;height:420px;"></object>
```

## Paraboloid
```@example geometricmodeling
plotly()
p = 2
k = KnotVector([-1,-1,-1,1,1,1])
P = BSplineSpace{p}(k)
a = [SVector(i,j,2i^2+2j^2-2) for i in -1:1, j in -1:1]
M = BSplineManifold(a,P,P)
plot(M)
savefig("geometricmodeling-paraboloid.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-paraboloid.html" style="width:100%;height:420px;"></object>
```

## Hyperbolic paraboloid
```@example geometricmodeling
plotly()
a = [SVector(i,j,2i^2-2j^2) for i in -1:1, j in -1:1]
M = BSplineManifold(a,P,P)
plot(M)
savefig("geometricmodeling-hyperbolicparaboloid.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../geometricmodeling-hyperbolicparaboloid.html" style="width:100%;height:420px;"></object>
```
