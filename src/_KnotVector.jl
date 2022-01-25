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

KnotVector(k::KnotVector) = k
KnotVector(k::AbstractKnotVector{T}) where T = unsafe_knotvector(T,_vec(k))
KnotVector{T}(k::KnotVector{T}) where T = k
KnotVector{T}(k::KnotVector{S}) where {T,S} = unsafe_knotvector(T,k.vector)
KnotVector{T}(k::AbstractKnotVector{S}) where {T,S} = unsafe_knotvector(T,_vec(k))

Base.convert(T::Type{<:KnotVector},k::AbstractKnotVector) = T(k)
function Base.promote_rule(::Type{KnotVector{T}}, ::Type{KnotVector{S}}) where {T,S}
    KnotVector{promote_type(T,S)}
end

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
KnotVector{T}() where T<:Real = unsafe_knotvector(T, T[])

function Base.show(io::IO, k::KnotVector)
    if k.vector == Float64[]
        print(io, "KnotVector([])")
    else
        print(io, "KnotVector($(k.vector))")
    end
end

"""
Convert `AbstractKnotVector` to `AbstractVector`
"""
_vec

_vec(k::KnotVector) = k.vector

Base.zero(::Type{<:KnotVector}) = KnotVector()
Base.zero(::KnotVector{T}) where T = KnotVector{T}()
Base.:(==)(k₁::AbstractKnotVector, k₂::AbstractKnotVector) = (_vec(k₁) == _vec(k₂))

Base.eltype(::AbstractKnotVector{T}) where T = T

# TODO: remove collect
Base.collect(k::KnotVector) = copy(k.vector)

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
Base.:+(k1::AbstractKnotVector, k2::AbstractKnotVector) = +(promote(k1,k2)...)

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
function Base.:*(m::Integer, k::AbstractKnotVector)
    if m == 0
        return zero(k)
    elseif m > 0
        return sum(k for _ in 1:m)
    else
        throw(DomainError(m, "The number to be multiplied must be non-negative."))
    end
end
Base.:*(k::AbstractKnotVector, m::Integer) = m*k

Base.in(r::Real, k::AbstractKnotVector) = in(r, _vec(k))
Base.getindex(k::AbstractKnotVector, i::Integer) = _vec(k)[i]
Base.getindex(k::KnotVector, v::AbstractVector{<:Integer}) = KnotVector(_vec(k)[v])

@doc raw"""
Length of knot vector

```math
\begin{aligned}
\#{k}
&=(\text{number of knot elements of} \  k) \\
\end{aligned}
```

For example, ``\#{(1,2,2,3)}=4``.

# Examples
```jldoctest
julia> k = KnotVector(1,2,2,3);

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
Check a inclusive relationship ``k\subseteq k'``, for example:
```math
\begin{aligned}
(1,2) &\subseteq (1,2,3) \\
(1,2,2) &\not\subseteq (1,2,3) \\
(1,2,3) &\subseteq (1,2,3) \\
\end{aligned}
```

# Examples
```jldoctest
julia> KnotVector(1,2) ⊆ KnotVector(1,2,3)
true

julia> KnotVector(1,2,2) ⊆ KnotVector(1,2,3)
false

julia> KnotVector(1,2,3) ⊆ KnotVector(1,2,3)
true
```
"""
function Base.issubset(k::KnotVector, k′::KnotVector)
    v = k′.vector
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
