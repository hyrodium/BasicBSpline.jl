# B-spline space

```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

Before defining B-spline space, we'll define polynomial space with degree ``p``.

!!! tip "Def.  Polynomial space"
    Polynomial space with degree ``p``.
    ```math
    \mathcal{P}[p]
    =\left\{f:\mathbb{R}\to\mathbb{R}\ ;\ t\mapsto a_0+a_1t^1+\cdots+a_pt^p \  \left| \ %
        a_i\in \mathbb{R}
        \right.
    \right\}
    ```
    This space ``\mathcal{P}[p]`` is a ``(p+1)``-dimensional linear space.

Note that ``\{t\mapsto t^i\}_{0 \le i \le p}`` is a basis of ``\mathcal{P}[p]``, and also the set of [Bernstein polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) ``\{B_{(i,p)}\}_i`` is a basis of ``\mathcal{P}[p]``.

```math
\begin{aligned}
B_{(i,p)}(t)
&=\binom{p}{i-1}t^{i-1}(1-t)^{p-i+1}
&(i=1, \dots, p+1)
\end{aligned}
```

Where ``\binom{p}{i-1}`` is a [binomial coefficient](https://en.wikipedia.org/wiki/Binomial_coefficient).

!!! tip "Def.  B-spline space"
    For given polynomial degree ``p\ge 0`` and knot vector ``k=(k_1,\dots,k_l)``, B-spline space ``\mathcal{P}[p,k]`` is defined as follows:
    ```math
    \mathcal{P}[p,k]
    =\left\{f:\mathbb{R}\to\mathbb{R} \  \left| \ %
        \begin{gathered}
            \operatorname{supp}(f)\subseteq [k_1, k_l] \\
            \exists \tilde{f}\in\mathcal{P}[p], f|_{[k_{i}, k_{i+1})} = \tilde{f}|_{[k_{i}, k_{i+1})}  \\
            \forall t \in \mathbb{R}, \exists \delta > 0, f|_{(t-\delta,t+\delta)}\in C^{p-\mathfrak{n}_k(t)}
        \end{gathered} \right.
    \right\}
    ```

Note that each element of the space ``\mathcal{P}[p,k]`` is a piecewise polynomial.

[TODO: fig]

```@repl math
p = 2
k = KnotVector([1,3,5,6,8,9])
BSplineSpace{p}(k)
```

```@docs
BSplineSpace
```

```@docs
UniformBSplineSpace
```

A B-spline space is said to be **non-degenerate** if its degree and knotvector satisfies following property:
```math
\begin{aligned}
k_{i}&<k_{i+p+1} & (1 \le i \le l-p-1)
\end{aligned}
```

```@docs
isnondegenerate
```

```@docs
isdegenerate(P::AbstractBSplineSpace)
```

The B-spline space is a linear space, and if a B-spline space is non-degenerate, its dimension is calculated by:
```math
\dim(\mathcal{P}[p,k])=\sharp k - p -1
```

```@repl math
dim(BSplineSpace{2}(KnotVector([1,3,5,6,8,9])))
```

```@docs
dim
```
