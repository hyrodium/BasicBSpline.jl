# B-spline basis function

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

```@setup math2
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
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

The next figure shows the plot of B-spline basis functions.
You can manipulate these plots on [desmos graphing calculator](https://www.desmos.com/calculator/ql6jqgdabs)!

![](img/bsplinebasis.png)

!!! info "Thm.  Basis of B-spline space"
    The set of functions ``\{B_{(i,p,k)}\}_i`` is a basis of B-spline space ``\mathcal{P}[p,k]``.


```@repl math
p = 2
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P = BSplineSpace{p}(k)
plot([t->bsplinebasis₊₀(P,i,t) for i in 1:dim(P)], 0, 10, ylims=(0,1), label=false)
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
p = 2
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P = BSplineSpace{p}(k)
plot([t->bsplinebasis₋₀(P,i,t) for i in 1:dim(P)], 0, 10, ylims=(0,1), label=false)
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

```@docs
bsplinesupport
```

## Partition of unity
!!! info "Thm.  Partition of unity"
    Let ``B_{(i,p,k)}`` be a B-spline basis function, then the following equation is satisfied.
    ```math
    \begin{aligned}
    \sum_{i}B_{(i,p,k)}(t) &= 1 & (k_{p+1} \le t < k_{l-p}) \\
    0 \le B_{(i,p,k)}(t) &\le 1
    \end{aligned}
    ```

```@repl math
p = 2
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis₊₀(P,i,t) for i in 1:dim(P)), 0, 10, ylims=(0,1.1), label=false)
savefig("sumofbsplineplot.html") # hide
```

```@raw html
<object type="text/html" data="../sumofbsplineplot.html" style="width:100%;height:420px;"></object>
```

To satisfy the partition of unity on whole interval ``[1,8]``, sometimes more knots will be inserted to the endpoints of the interval.

```@repl math
p = 2
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]) + p * KnotVector([0,10])
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis₊₀(P,i,t) for i in 1:dim(P)), 0, 10, ylims=(0,1.1), label=false)
savefig("sumofbsplineplot2.html") # hide
```

```@raw html
<object type="text/html" data="../sumofbsplineplot2.html" style="width:100%;height:420px;"></object>
```

But, the sum ``\sum_{i} B_{(i,p,k)}(t)`` is not equal to ``1`` at ``t=8``.
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
p = 2
k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]) + p * KnotVector([0,10])
P = BSplineSpace{p}(k)
plot(t->sum(bsplinebasis(P,i,t) for i in 1:dim(P)), 0, 10, ylims=(0,1.1), label=false)
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

## B-spline basis functions at specific point
Sometimes, you may need the non-zero values of B-spline basis functions at specific point.
The `bsplinebasisall` function is much more efficient than evaluating B-spline functions one by one with `bsplinebasis` function.

```@repl
using BenchmarkTools, BasicBSpline
P = BSplineSpace{3}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))
t = 6.3
bsplinebasis.(P, 1:4, t)
bsplinebasisall(P, 1, t)
@benchmark bsplinebasis.($P, 1:4, $t)
@benchmark bsplinebasisall($P, 1, $t)
```

```@docs
intervalindex
```

```@docs
bsplinebasisall
```

The next figures illustlates the relation between `domain(P)`, `intervalindex(P,t)` and `bsplinebasisall(P,i,t)`.

```@example math2
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

## Uniform B-spline basis and uniform distribution

Let ``X_1, \dots, X_n`` be i.i.d. random variables with ``X_i \sim U(0,1)``, then the probability density function of ``X_1+\cdots+X_n`` can be obtained via `BasicBSpline.uniform_bsplinebasis_kernel(Val(n-1),t)`.

```@example math
N = 100000
# polynomial degree 0
plot1 = histogram([rand() for _ in 1:N], normalize=true, label=false)
plot!(t->BasicBSpline.uniform_bsplinebasis_kernel(Val(0),t), label=false)
# polynomial degree 1
plot2 = histogram([rand()+rand() for _ in 1:N], normalize=true, label=false)
plot!(t->BasicBSpline.uniform_bsplinebasis_kernel(Val(1),t), label=false)
# polynomial degree 2
plot3 = histogram([rand()+rand()+rand() for _ in 1:N], normalize=true, label=false)
plot!(t->BasicBSpline.uniform_bsplinebasis_kernel(Val(2),t), label=false)
# polynomial degree 3
plot4 = histogram([rand()+rand()+rand()+rand() for _ in 1:N], normalize=true, label=false)
plot!(t->BasicBSpline.uniform_bsplinebasis_kernel(Val(3),t), label=false)
# plot all
plot(plot1,plot2,plot3,plot4)
savefig("histogram-uniform.html") # hide
```

```@raw html
<object type="text/html" data="../histogram-uniform.html" style="width:100%;height:420px;"></object>
```
