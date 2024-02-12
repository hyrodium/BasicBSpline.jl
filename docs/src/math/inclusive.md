# Inclusive relation between B-spline spaces

## Setup

```@example math_inclusive
using BasicBSpline
using StaticArrays
using Plots
```

## Theorem on [`issubset`](@ref)

!!! info "Thm.  Inclusive relation between B-spline spaces"
    For non-degenerate B-spline spaces, the following relationship holds.
    ```math
    \mathcal{P}[p,k]
    \subseteq \mathcal{P}[p',k']
    \Leftrightarrow (m=p'-p \ge 0 \ \text{and} \ k+m\widehat{k}\subseteq k')
    ```

### Examples

Here are plots of the B-spline basis functions of the spaces `P1`, `P2`, `P3`.

```@example math_inclusive
P1 = BSplineSpace{1}(KnotVector([1,3,6,6]))
P2 = BSplineSpace{1}(KnotVector([1,3,5,6,6,8,9]))
P3 = BSplineSpace{2}(KnotVector([1,1,3,3,6,6,6,8,9]))
plotly()
plot(
    plot(P1; ylims=(0,1), label="P1"),
    plot(P2; ylims=(0,1), label="P2"),
    plot(P3; ylims=(0,1), label="P3"),
    layout=(3,1),
    link=:x
)
savefig("inclusive-issubset.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../inclusive-issubset.html" style="width:100%;height:420px;"></object>
```

These spaces have the folllowing incusive relationships.

```@repl math_inclusive
P1 ⊆ P2
P1 ⊆ P3
P2 ⊆ P3, P3 ⊆ P2, P2 ⊆ P1, P3 ⊆ P1
```

## Definition on [`issqsubset`](@ref)

!!! tip "Def.  Inclusive relation between B-spline spaces"
    For non-degenerate B-spline spaces, the following relationship holds.
    ```math
    \mathcal{P}[p,k]
    \sqsubseteq\mathcal{P}[p',k']
    \Leftrightarrow
    \mathcal{P}[p,k]|_{[k_{p+1},k_{l-p}]}
    \subseteq\mathcal{P}[p',k']|_{[k'_{p'+1},k'_{l'-p'}]}
    ```

### Examples

Here are plots of the B-spline basis functions of the spaces `P1`, `P2`, `P3`.

```@example math_inclusive
P1 = BSplineSpace{1}(KnotVector([1,3,6,6]))  # Save definition as above
P4 = BSplineSpace{1}(KnotVector([1,3,5,6,7]))
P5 = BSplineSpace{2}(KnotVector([1,1,3,3,6,6,7,9]))
plotly()
plot(
    plot(P1; ylims=(0,1), legend=false),
    plot(P4; ylims=(0,1), legend=false),
    plot(P5; ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)
savefig("inclusive-issqsubset.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../inclusive-issqsubset.html" style="width:100%;height:420px;"></object>
```

These spaces have the folllowing incusive relationships.

```@repl math_inclusive
P1 ⊑ P4
P1 ⊑ P5
P4 ⊑ P5, P5 ⊑ P4, P4 ⊑ P1, P5 ⊑ P1
```

## Change basis with a matrix

If ``\mathcal{P}[p,k] \subseteq \mathcal{P}[p',k']`` (or ``\mathcal{P}[p,k] \sqsubseteq \mathcal{P}[p',k']``), there exists a ``n \times n'`` matrix ``A`` which holds:

```math
\begin{aligned}
B_{(i,p,k)}
&=\sum_{j}A_{ij} B_{(j,p',k')} \\
n&=\dim(\mathcal{P}[p,k]) \\
n'&=\dim(\mathcal{P}[p',k'])
\end{aligned}
```

You can calculate the change of basis matrix ``A`` with [`changebasis`](@ref).

```@repl math_inclusive
A12 = changebasis(P1,P2)
A13 = changebasis(P1,P3)
```

```@example math_inclusive
plot(
    plot([t->bsplinebasis₊₀(P1,i,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    plot([t->sum(A12[i,j]*bsplinebasis₊₀(P2,j,t) for j in 1:dim(P2)) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    plot([t->sum(A13[i,j]*bsplinebasis₊₀(P3,j,t) for j in 1:dim(P3)) for i in 1:dim(P1)], 1, 9, ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)
savefig("inclusive-issubset-matrix.html") # hide
nothing # hide
```

```@raw html
<object type="text/html" data="../inclusive-issubset-matrix.html" style="width:100%;height:420px;"></object>
```

## Expand spaces with additional knots or polynomial degree

There are some functions to expand spaces with additional knots or polynomial degree.

* [`expandspace_R`](@ref)
* [`expandspace_I`](@ref)
* [`expandspace`](@ref)

```@repl math_inclusive
P = BSplineSpace{2}(knotvector"21 113")
P_R = expandspace_R(P, Val(1), KnotVector([3.4, 4.2]))
P_I = expandspace_I(P, Val(1), KnotVector([3.4, 4.2]))
P ⊆ P_R
P ⊆ P_I
P ⊑ P_R
P ⊑ P_I
```
