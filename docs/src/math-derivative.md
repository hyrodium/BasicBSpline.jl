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

```@example math
using Plots; plotlyjs()
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P = BSplineSpace{3}(k)
plot(
    plot(BSplineDerivativeSpace{0}(P), label="0th derivative", color=:black),
    plot(BSplineDerivativeSpace{1}(P), label="1st derivative", color=:red),
    plot(BSplineDerivativeSpace{2}(P), label="2nd derivative", color=:green),
    plot(BSplineDerivativeSpace{3}(P), label="3rd derivative", color=:blue),
)
savefig("bsplinebasisderivativeplot.html") # hide
```

```@raw html
<object type="text/html" data="../bsplinebasisderivativeplot.html" style="width:100%;height:420px;"></object>
```

```@docs
BSplineDerivativeSpace
```

```@docs
bsplinebasis′₊₀
```

```@docs
bsplinebasis′₋₀
```

```@docs
bsplinebasis′
```
