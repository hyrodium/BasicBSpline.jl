# Inclusive relation between B-spline spaces

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

!!! info "Thm.  Inclusive relation between B-spline spaces"
    For non-degenerate B-spline spaces, the following relationship holds.
    ```math
    \mathcal{P}[p,k]
    \subseteq \mathcal{P}[p',k']
    \Leftrightarrow (m=p'-p \ge 0 \ \text{and} \ k+m\widehat{k}\subseteq k')
    ```

```@docs
Base.issubset(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
```

Here are plots of the B-spline basis functions of the spaces `P1`, `P2`, `P3`.

```@repl math
P1 = BSplineSpace{1}(KnotVector([1,3,5,8]))
P2 = BSplineSpace{1}(KnotVector([1,3,5,6,8,9]))
P3 = BSplineSpace{2}(KnotVector([1,1,3,3,5,5,8,8]))
plot(
    plot([t->bsplinebasis₊₀(P1,i,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(P2,i,t) for i in 1:dim(P2)], 1, 9, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(P3,i,t) for i in 1:dim(P3)], 1, 9, ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)
savefig("subbsplineplot.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../subbsplineplot.html" style="width:100%;height:420px;"></object>
```

This means, there exists a ``n \times n'`` matrix ``A`` which holds:

```math
\begin{aligned}
B_{(i,p,k)}
&=\sum_{j}A_{ij} B_{(j,p',k')} \\
n&=\dim(\mathcal{P}[p,k]) \\
n'&=\dim(\mathcal{P}[p',k'])
\end{aligned}
```

You can calculate the change of basis matrix ``A`` with `changebasis`.

```@repl math
A12 = changebasis(P1,P2)
A13 = changebasis(P1,P3)
```

```@repl math
plot(
    plot([t->bsplinebasis₊₀(P1,i,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    plot([t->sum(A12[i,j]*bsplinebasis₊₀(P2,j,t) for j in 1:dim(P2)) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    plot([t->sum(A13[i,j]*bsplinebasis₊₀(P3,j,t) for j in 1:dim(P3)) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)
savefig("subbsplineplot2.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../subbsplineplot2.html" style="width:100%;height:420px;"></object>
```

```@docs
issqsubset
```

```@docs
expandspace
```

```@docs
expandspace_R
```

```@docs
expandspace_I
```
