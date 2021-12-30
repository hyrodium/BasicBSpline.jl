# KnotVector

abstract type AbstractKnotVector{T<:Real} end

@doc raw"""
Construct knot vector from given array.
```math
k=(k_1,\dots,k_l)
```

# Examples
```jldoctest
julia> k = KnotVector([1,2,3])
KnotVector([1.0, 2.0, 3.0])

julia> k = KnotVector(1:3)
KnotVector([1.0, 2.0, 3.0])
```
"""
struct KnotVector{T} <: AbstractKnotVector{T}
    vector::Vector{T}
    global unsafe_knotvector(::Type{T}, v) where T = new{T}(v)
end
KnotVector{T}(v::AbstractVector) where T = unsafe_knotvector(T,sort(v))
KnotVector(v::AbstractVector{T}) where {T<:Real} = unsafe_knotvector(float(T),sort(v))

@doc raw"""
Construct knot vector from given real numbers.

# Examples
```jldoctest
julia> k = KnotVector(1,2,3)
KnotVector([1.0, 2.0, 3.0])

julia> k = KnotVector()
KnotVector([])
```
"""
function KnotVector(knots::Real...)
    return KnotVector(collect(knots))
end
function KnotVector{T}(knots::Real...) where T<:Real
    return unsafe_knotvector(T, sort!(collect(knots)))
end
KnotVector() = unsafe_knotvector(Float64, Float64[])

function Base.show(io::IO, k::KnotVector)
    if k.vector == Float64[]
        print(io, "KnotVector([])")
    else
        print(io, "KnotVector($(k.vector))")
    end
end

Base.zero(::Type{<:KnotVector}) = KnotVector()
Base.:(==)(k₁::KnotVector, k₂::KnotVector) = (k₁.vector == k₂.vector)

@doc raw"""
Sum of knot vectors

```math
\begin{aligned}
k^{(1)}+k^{(2)}
&=(k^{(1)}_1, \dots, k^{(1)}_{l^{(1)}}) + (k^{(2)}_1, \dots, k^{(2)}_{l^{(2)}}) \\
&=(\text{sort of union of} \  k^{(1)} \ \text{and} \  k^{(2)} \text{)}
\end{aligned}
```

For example, ``(1,2,3,5)+(4,5,8)=(1,2,3,4,5,5,8)``.

# Examples
```jldoctest
julia> k1 = KnotVector(1,2,3,5);

julia> k2 = KnotVector(4,5,8);

julia> k1 + k2
KnotVector([1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 8.0])
```
"""
Base.:+(k1::KnotVector{T}, k2::KnotVector{T}) where T = unsafe_knotvector(T,sort!(vcat(k1.vector,k2.vector)))

# TODO: add a method for KnotVector{Int} + KnotVector{Float64}

@doc raw"""
Product of integer and knot vector

```math
\begin{aligned}
m\cdot k&=\underbrace{k+\cdots+k}_{m}
\end{aligned}
```

For example, ``2\cdot (1,2,2,5)=(1,1,2,2,2,2,5,5)``.

# Examples
```jldoctest
julia> k = KnotVector(1,2,2,5);

julia> 2 * k
KnotVector([1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0])
```
"""
function Base.:*(p₊::Integer, k::AbstractKnotVector)
    if p₊ == 0
        return zero(KnotVector)
    elseif p₊ > 0
        return sum(k for _ in 1:p₊)
    else
        throw(DomainError(p₊, "additional degree of polynominal must be non-negative"))
    end
end
Base.:*(k::AbstractKnotVector, p₊::Integer) = p₊*k

Base.in(r::Real, k::KnotVector) = in(r, k.vector)
Base.getindex(k::KnotVector, i::Integer) = k.vector[i]
Base.getindex(k::KnotVector, v::AbstractVector{<:Integer}) = KnotVector(k.vector[v])

@doc raw"""
Length of knot vector

# Examples
```jldoctest
julia> k = KnotVector(1,2,3,5);

julia> length(k)
4
```
"""
Base.length(k::KnotVector) = length(k.vector)

Base.firstindex(k::KnotVector) = 1
Base.lastindex(k::KnotVector) = length(k)

@doc raw"""
Unique elements of knot vector.

```math
\begin{aligned}
\widehat{k}
&=(\text{unique knot elements of} \  k) \\
\end{aligned}
```

For example, ``\widehat{(1,2,2,3)}=(1,2,3)``.

# Examples
```jldoctest
julia> k = KnotVector([1,2,2,3]);

julia> unique(k)
KnotVector([1.0, 2.0, 3.0])
```
"""
Base.unique(k::KnotVector) = KnotVector(unique(k.vector))
Base.unique!(k::KnotVector) = KnotVector(unique!(k.vector))
Base.iterate(k::KnotVector) = iterate(k.vector)
Base.iterate(k::KnotVector, i::Integer) = iterate(k.vector, i)
Base.searchsortedfirst(k::KnotVector,t) = searchsortedfirst(k.vector,t)
Base.searchsortedlast(k::KnotVector,t) = searchsortedlast(k.vector,t)
Base.searchsorted(k::KnotVector,t) = searchsorted(k.vector,t)

@doc raw"""
Check a inclusive relation ship ``k\subset k'``.
```math
(1,2,3) \subseteq (1,2,3,4)
(1,2,5) \nsubseteq (1,2,3,4)
```
"""
function Base.issubset(k::KnotVector, k′::KnotVector)
    v = k′.vector
    l = length(v)
    i = 0
    for kᵢ in k
        i = findnext(==(kᵢ), v, i+1)
        if isnothing(i)
            return false
        end
    end
    return true
end

@doc raw"""
For given knot vector ``k``, the following function ``\mathfrak{n}_k:\mathbb{R}\to\mathbb{Z}`` represents the number of knots that duplicate the knot vector ``k``.

```math
\mathfrak{n}_k(t) = \sharp\{i \mid k_i=t \}
```
For example, if ``k=(1,2,2,3)``, then ``\mathfrak{n}_k(0.3)=0``, ``\mathfrak{n}_k(1)=1``, ``\mathfrak{n}_k(2)=2``.

```jldoctest
julia> k = KnotVector([1,2,2,3]);

julia> 𝔫(k,0.3)
0

julia> 𝔫(k,1.0)
1

julia> 𝔫(k,2.0)
2
```
"""
function 𝔫(k::KnotVector, t::Real)
    # for small case, this is faster
    # return count(==(t), k.vector)

    # for large case, this is faster
    return length(searchsorted(k,t))
end
