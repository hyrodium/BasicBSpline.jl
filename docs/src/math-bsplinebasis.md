# B-spline basis function

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots
```

## Basic properties of B-spline basis function

!!! tip "Def.  B-spline space"
    B-spline basis function is defined by [Cox–de Boor recursion formula](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm).
    ```math
    \begin{aligned}
    {B}_{(i,p,k)}(t)
    &=
    \frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
    +\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
    {B}_{(i,0,k)}(t)
    &=
    \begin{cases}
        &1\quad (k_{i}\le t< k_{i+1})\\
        &0\quad (\text{otherwise})
    \end{cases}
    \end{aligned}
    ```
    If the denominator is zero, then the term is assumed to be zero.


!!! info "Thm.  Basis of B-spline space"
    The set of functions ``\{B_{(i,p,k)}\}_i`` is a basis of B-spline space ``\mathcal{P}[p,k]``.


```@repl math
using Plots; plotly()
p = 2
k = KnotVector(1:8)
P = BSplineSpace{p}(k)
plot([t->bsplinebasis₊₀(P,i,t) for i in 1:dim(P)], 1, 8, ylims=(0,1), label=false)
savefig("bsplinebasisplot.html") # hide
```

```@raw html
<object type="text/html" data="../bsplinebasisplot.html" style="width:100%;height:420px;"></object>
```

You can choose the first terms in different ways.

```math
\begin{aligned}
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i} < t \le k_{i+1}) \\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

```@repl math
using Plots; plotly()
p = 2
k = KnotVector(1:8)
P = BSplineSpace{p}(k)
plot([t->bsplinebasis₋₀(P,i,t) for i in 1:dim(P)], 1, 8, ylims=(0,1), label=false)
savefig("bsplinebasisplot2.html") # hide
```

```@raw html
<object type="text/html" data="../bsplinebasisplot2.html" style="width:100%;height:420px;"></object>
```

In these cases, each B-spline basis function ``B_{(i,2,k)}`` is coninuous, so `bsplinebasis₊₀` and `bsplinebasis₋₀` are equal.

## Support of B-spline basis function
!!! info "Thm.  Support of B-spline basis function"
    If a B-spline space``\mathcal{P}[p,k]`` is non-degenerate, the support of its basis function is calculated as follows:
    ```math
    \operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
    ```

[TODO: fig]

```@repl math
i = 2
k = KnotVector([5,12,13,13,14])
p = 2
P = BSplineSpace{p}(k)
bsplinesupport(P,i) # 12..14
```

## Partition of unity
!!! info "Thm.  Partition of unity"
    ```math
    \begin{aligned}
    \sum_{i}B_{(i,p,k)}(t) &= 1 & (k_{p+1} \le t < k_{l-p}) \\
    0 \le B_{(i,p,k)}(t) &\le 1
    \end{aligned}
    ```

```@repl math
using Plots; plotly()
p = 2
k = KnotVector(1:8)
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis₊₀(P,i,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.1), label=false)
savefig("sumofbsplineplot.html") # hide
```

```@raw html
<object type="text/html" data="../sumofbsplineplot.html" style="width:100%;height:420px;"></object>
```

To satisfy the partition of unity on whole interval ``[1,8]``, sometimes more knots will be inserted to the endpoints of the interval.

```@repl math
using Plots; plotly()
p = 2
k = KnotVector(1:8) + p * KnotVector([1,8])
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis₊₀(P,i,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.1), label=false)
savefig("sumofbsplineplot2.html") # hide
```

```@raw html
<object type="text/html" data="../sumofbsplineplot2.html" style="width:100%;height:420px;"></object>
```

But, the sum ``\sum_{i} B_{(i,p,k)}(t)`` is not equal to ``1`` if ``t=8``.
Therefore, to satisfy partition of unity on closed interval ``[k_{p+1}, k_{l-p}]``, the definition of first terms of B-spline basis functions are sometimes replaced:

```math
\begin{aligned}
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i} \le t<k_{i+1})\\
    &1\quad (k_{i} < t = k_{i+1}=k_{l})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

```@repl math
using Plots; plotly()
p = 2
k = KnotVector(1:8) + p * KnotVector([1,8])
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis(P,i,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.1), label=false)
savefig("sumofbsplineplot3.html") # hide
```

```@raw html
<object type="text/html" data="../sumofbsplineplot3.html" style="width:100%;height:420px;"></object>
```

```@docs
bsplinebasis₊₀
```

```@docs
bsplinebasis₋₀
```

```@docs
bsplinebasis
```

```@docs
bsplinesupport
```

```@docs
intervalindex
```

```@docs
bsplinebasisall
```

The next figures illustlates the relation between `domain(P)`, `intervalindex(P,t)` and `bsplinebasisall(P,i,t)`.

```@example
using BasicBSpline
using Plots; plotly()

k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

for p in 1:3
    P = BSplineSpace{p}(k)
    plot(P, legend=:topleft, label="B-spline basis (p=1)")
    plot!(t->intervalindex(P,t),0,10, label="Interval index")
    plot!(t->sum(bsplinebasis(P,i,t) for i in 1:dim(P)),0,10, label="Sum of B-spline basis")
    scatter!(k.vector,zero(k.vector), label="knot vector")
    plot!([t->bsplinebasisall(P,1,t)[i] for i in 1:p+1],0,10, color=:black, label="bsplinebasisall (i=1)", ylim=(-1,8-2p))
    savefig("bsplinebasisall-$(p).html") # hide
end
```

```@raw html
<object type="text/html" data="../bsplinebasisall-1.html" style="width:100%;height:420px;"></object>
```

```@raw html
<object type="text/html" data="../bsplinebasisall-2.html" style="width:100%;height:420px;"></object>
```

```@raw html
<object type="text/html" data="../bsplinebasisall-3.html" style="width:100%;height:420px;"></object>
```
