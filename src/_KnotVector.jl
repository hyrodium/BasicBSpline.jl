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
KnotVector([1, 2, 3])

julia> k = KnotVector(1:3)
KnotVector([1, 2, 3])
```
"""
struct KnotVector{T} <: AbstractKnotVector{T}
    vector::Vector{T}
    global unsafe_knotvector(::Type{T}, v) where T = new{T}(v)
end
KnotVector{T}(v::AbstractVector) where T = unsafe_knotvector(T,sort(v))
KnotVector(v::AbstractVector{T}) where T = unsafe_knotvector(T,sort(v))
AbstractKnotVector{S}(k::KnotVector{T}) where {S, T} = unsafe_knotvector(promote_type(S,T), _vec(k))

Base.copy(k::KnotVector{T}) where T = unsafe_knotvector(T,copy(_vec(k)))

KnotVector(k::KnotVector) = k
KnotVector(k::AbstractKnotVector{T}) where T = unsafe_knotvector(T,_vec(k))
KnotVector{T}(k::KnotVector{T}) where T = k
KnotVector{T}(k::AbstractKnotVector) where T = unsafe_knotvector(T,_vec(k))

Base.convert(T::Type{<:AbstractKnotVector}, k::AbstractKnotVector) = T(k)
function Base.promote_rule(::Type{KnotVector{T}}, ::Type{KnotVector{S}}) where {T,S}
    KnotVector{promote_type(T,S)}
end

function Base.show(io::IO, k::KnotVector)
    print(io, "KnotVector($(k.vector))")
end

function Base.show(io::IO, k::T) where T<:AbstractKnotVector
    print(io, "$(nameof(T))($(k.vector))")
end

"""
Convert `AbstractKnotVector` to `AbstractVector`
"""
_vec

_vec(k::KnotVector) = k.vector

Base.zero(::Type{KnotVector}) = zero(KnotVector{Float64})
Base.zero(::Type{KnotVector{T}}) where T = unsafe_knotvector(T, T[])
Base.zero(::KnotVector{T}) where T = unsafe_knotvector(T, T[])
Base.:(==)(k₁::AbstractKnotVector, k₂::AbstractKnotVector) = (_vec(k₁) == _vec(k₂))

Base.eltype(::AbstractKnotVector{T}) where T = T

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
julia> k1 = KnotVector([1,2,3,5]);

julia> k2 = KnotVector([4,5,8]);

julia> k1 + k2
KnotVector([1, 2, 3, 4, 5, 5, 8])
```
"""
function Base.:+(k1::KnotVector{T}, k2::KnotVector{T}) where T
    v1, v2 = k1.vector, k2.vector
    n1, n2 = length(v1), length(v2)
    iszero(n1) && return k2
    iszero(n2) && return k1
    n = n1 + n2
    v = Vector{T}(undef, n)
    i1 = i2 = 1
    for i in 1:n
        if isless(v1[i1], v2[i2])
            v[i] = v1[i1]
            i1 += 1
            if i1 > n1
                v[i+1:n] = view(v2, i2:n2)
                break
            end
        else
            v[i] = v2[i2]
            i2 += 1
            if i2 > n2
                v[i+1:n] = view(v1, i1:n1)
                break
            end
        end
    end
    return BasicBSpline.unsafe_knotvector(T,v)
end
Base.:+(k1::AbstractKnotVector{T1}, k2::AbstractKnotVector{T2}) where {T1, T2} = +(promote(k1,k2)...)

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
julia> k = KnotVector([1,2,2,5]);

julia> 2 * k
KnotVector([1, 1, 2, 2, 2, 2, 5, 5])
```
"""
function Base.:*(m::Integer, k::AbstractKnotVector)
    return m*KnotVector(k)
end
function Base.:*(m::Integer, k::KnotVector{T}) where T
    if m == 0
        return zero(k)
    elseif m == 1
        return k
    elseif m > 1
        n = length(k)
        v = Vector{T}(undef, m*n)
        for i in 1:m
            v[i:m:m*(n-1)+i] .= k.vector
        end
        return BasicBSpline.unsafe_knotvector(T,v)
    else
        throw(DomainError(m, "The number to be multiplied must be non-negative."))
    end
end
Base.:*(k::AbstractKnotVector, m::Integer) = m*k

Base.in(r::Real, k::AbstractKnotVector) = in(r, _vec(k))
Base.getindex(k::AbstractKnotVector, i::Integer) = _vec(k)[i]
Base.getindex(k::AbstractKnotVector, v::AbstractVector{<:Integer}) = KnotVector(_vec(k)[v])

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
julia> k = KnotVector([1,2,2,3]);

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
KnotVector([1, 2, 3])
```
"""
Base.unique(k::KnotVector) = KnotVector(unique(k.vector))
Base.unique!(k::KnotVector) = KnotVector(unique!(k.vector))
Base.iterate(k::AbstractKnotVector) = iterate(_vec(k))
Base.iterate(k::AbstractKnotVector, i) = iterate(_vec(k), i)
Base.searchsortedfirst(k::AbstractKnotVector,t) = searchsortedfirst(_vec(k),t)
Base.searchsortedlast(k::AbstractKnotVector,t) = searchsortedlast(_vec(k),t)
Base.searchsorted(k::AbstractKnotVector,t) = searchsorted(_vec(k),t)

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
julia> KnotVector([1,2]) ⊆ KnotVector([1,2,3])
true

julia> KnotVector([1,2,2]) ⊆ KnotVector([1,2,3])
false

julia> KnotVector([1,2,3]) ⊆ KnotVector([1,2,3])
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

function Base.issubset(k::AbstractKnotVector, k′::AbstractKnotVector)
    issubset(KnotVector(k),KnotVector(k′))
end

Base.:⊊(A::AbstractKnotVector, B::AbstractKnotVector) = (A ≠ B) & (A ⊆ B)
Base.:⊋(A::AbstractKnotVector, B::AbstractKnotVector) = (A ≠ B) & (A ⊇ B)

@doc raw"""
For given knot vector ``k``, the following function ``\mathfrak{n}_k:\mathbb{R}\to\mathbb{Z}`` represents the number of knots that duplicate the knot vector ``k``.

```math
\mathfrak{n}_k(t) = \#\{i \mid k_i=t \}
```
For example, if ``k=(1,2,2,3)``, then ``\mathfrak{n}_k(0.3)=0``, ``\mathfrak{n}_k(1)=1``, ``\mathfrak{n}_k(2)=2``.

```jldoctest
julia> k = KnotVector([1,2,2,3]);

julia> countknots(k,0.3)
0

julia> countknots(k,1.0)
1

julia> countknots(k,2.0)
2
```
"""
function countknots(k::AbstractKnotVector, t::Real)
    # for small case, this is faster
    # return count(==(t), k.vector)

    # for large case, this is faster
    return length(searchsorted(_vec(k), t, lt=<))
end

function Base.float(k::KnotVector{T}) where T
    KnotVector(float(_vec(k)))
end
