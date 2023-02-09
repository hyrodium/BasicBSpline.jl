# Knot vector

## Definition
```@setup math
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using Plots; plotly()
```

!!! tip "Def.  Knot vector"
    A finite sequence
    ```math
    k = (k_1, \dots, k_l)
    ```
    is called **knot vector** if the sequence is broad monotonic increase, i.e. ``k_{i} \le k_{i+1}``.


[TODO: fig]

```@docs
KnotVector
```

```@docs
UniformKnotVector
```

```@docs
EmptyKnotVector
```

## Operations for knot vectors

```@docs
length(k::AbstractKnotVector)
```

Although a knot vector is **not** a vector in linear algebra, but we introduce **additional operator** ``+``.

```@docs
Base.:+(k1::KnotVector{T}, k2::KnotVector{T}) where T
```

Note that the operator `+(::KnotVector, ::KnotVector)` is commutative.
This is why we choose the ``+`` sign.
We also introduce **product operator** ``\cdot`` for knot vector.

```@docs
*(m::Integer, k::AbstractKnotVector)
```

Inclusive relationship between knot vectors.

```@docs
Base.issubset(k::KnotVector, k′::KnotVector)
```

```@docs
unique(k::AbstractKnotVector)
```

```@docs
countknots(k::AbstractKnotVector, t::Real)
```
