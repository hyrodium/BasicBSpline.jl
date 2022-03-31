# Derivative of B-spline

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots
```

!!! info "Thm.  Derivative of B-spline basis function"
    The derivative of B-spline basis function can be expressed as follows:
    ```math
    \begin{aligned}
    \dot{B}_{(i,p,k)}(t)
    &=\frac{d}{dt}B_{(i,p,k)}(t) \\
    &=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
    \end{aligned}
    ```
    Note that ``\dot{B}_{(i,p,k)}\in\mathcal{P}[p-1,k]``.

```@repl math
using Plots
plotlyjs() # hide
p = 2
k = KnotVector(1:8)
P = BSplineSpace{p}(k)
plot([t->bsplinebasis′₊₀(P,i,t) for i in 1:dim(P)], 1, 8, ylims=(-2,2), label=false)
savefig("bsplinebasisderivativeplot.html") # hide
```

```@raw html
<object type="text/html" data="../bsplinebasisderivativeplot.html" style="width:100%;height:420px;"></object>
```
