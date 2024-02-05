# Knot vector

## Setup
```@example math_knotvector
using BasicBSpline
using InteractiveUtils  # hide
using Plots; gr()
```

## Definition

!!! tip "Def.  Knot vector"
    A finite sequence
    ```math
    k = (k_1, \dots, k_l)
    ```
    is called **knot vector** if the sequence is broad monotonic increase, i.e. ``k_{i} \le k_{i+1}``.

There are four sutypes of [`AbstractKnotVector`](@ref); [`KnotVector`](@ref), [`UniformKnotVector`](@ref), [`SubKnotVector`](@ref), and [`EmptyKnotVector`](@ref).

```@repl math_knotvector
subtypes(AbstractKnotVector)
```

They can be constructed like this.

```@repl math_knotvector
KnotVector([1,2,3])  # `KnotVector` stores a vector with `Vector{<:Real}`
KnotVector(1:3)
UniformKnotVector(1:8)  # `UniformKnotVector` stores a vector with `<:AbstractRange`
UniformKnotVector(8:-1:3)
view(KnotVector([1,2,3]), 2:3)  # Efficient slicing with `view`
EmptyKnotVector()  # Sometimes `EmptyKnotVector` is useful.
EmptyKnotVector{Float64}()
```

There is a useful string macro [`@knotvector_str`](@ref) that generates a `KnotVector` instance.

```@repl math_knotvector
knotvector"1231123"
```

A knot vector can be visualized with `Plots.plot`.

```@repl math_knotvector
k1 = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
plot(k1; label="k1")
k2 = knotvector"1231123"
plot!(k2; label="k2", shift_y=0.5, ylims=(-0.1,1), xticks=0:10)
savefig("knotvector.png") # hide
nothing # hide
```

![](knotvector.png)

## Operations for knot vectors

### Length of a knot vector

```@repl math_knotvector
k = KnotVector([1,2,2,3]);
length(k)
```

### Addition of knot vectors

Although a knot vector is **not** a vector in linear algebra, but we introduce **additional operator** ``+``.

```@repl math_knotvector
k1 = KnotVector([1,2,3,5]);
k2 = KnotVector([4,5,8]);
k1 + k2
```

Note that the operator `+(::KnotVector, ::KnotVector)` is commutative.
This is why we choose the ``+`` sign.
We also introduce **product operator** ``\cdot`` for knot vector.

### Multiplication of knot vectors

```@repl math_knotvector
k = KnotVector([1,2,2,5]);
2 * k
```

### Inclusive relationship between knot vectors

```@repl math_knotvector
KnotVector([1,2]) ⊆ KnotVector([1,2,3])
KnotVector([1,2,2]) ⊆ KnotVector([1,2,3])
KnotVector([1,2,3]) ⊆ KnotVector([1,2,3])
```

### Generate a knot vector with unique elements

```@repl math_knotvector
k = KnotVector([1,2,2,3]);
unique(k)
```

### Count knots in a knot vector

```@repl math_knotvector
k = KnotVector([1,2,2,3]);
countknots(k,0.3)
countknots(k,1.0)
countknots(k,2.0)
```
