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
view(KnotVector([1,2,3]), 2:3)  # Efficient and lazy knot vector with `view`
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
plot!(k2; label="k2", offset=0.5, ylims=(-0.1,1), xticks=0:10)
savefig("knotvector.png") # hide
nothing # hide
```

![](knotvector.png)

## Operations for knot vectors

### Setup and visualization

```@example math_knotvector
k1 = knotvector"1 2 1 11 2 31"
k2 = knotvector"113 1 12 2 32 1  1"
plot(k1; offset=-0.0, label="k1", xticks=1:18, yticks=nothing, legend=:right)
plot!(k2; offset=-0.2, label="k2")
plot!(k1+k2; offset=-0.4, label="k1+k2")
plot!(2k1; offset=-0.6, label="2k1")
plot!(unique(k1); offset=-0.8, label="unique(k1)")
plot!(unique(k2); offset=-1.0, label="unique(k2)")
plot!(unique(k1+k2); offset=-1.2, label="unique(k1+k2)")
savefig("knotvector_operations.png") # hide
nothing # hide
```

![](knotvector_operations.png)

### Length of a knot vector

[`length(k::AbstractKnotVector)`](@ref)

```@repl math_knotvector
length(k1)
length(k2)
```

### Addition of knot vectors

Although a knot vector is **not** a vector in linear algebra, but we introduce **additional operator** ``+``.

[`Base.:+(k1::KnotVector{T}, k2::KnotVector{T}) where T`](@ref)

```@repl math_knotvector
k1 + k2
```

Note that the operator `+(::KnotVector, ::KnotVector)` is commutative.
This is why we choose the ``+`` sign.
We also introduce **product operator** ``\cdot`` for knot vector.

### Multiplication of knot vectors

[`*(m::Integer, k::AbstractKnotVector)`](@ref)

```@repl math_knotvector
2*k1
2*k2
```

### Generate a knot vector with unique elements

[`unique(k::AbstractKnotVector)`](@ref)

```@repl math_knotvector
unique(k1)
unique(k2)
```

### Inclusive relationship between knot vectors

[`Base.issubset(k::KnotVector, k′::KnotVector)`](@ref)

```@repl math_knotvector
unique(k1) ⊆ k1 ⊆ k2
k1 ⊆ k1
k2 ⊆ k1
```

### Count knots in a knot vector

[`countknots(k::AbstractKnotVector, t::Real)`](@ref)

```@repl math_knotvector
countknots(k1, 0.5)
countknots(k1, 1.0)
countknots(k1, 3.0)
```
