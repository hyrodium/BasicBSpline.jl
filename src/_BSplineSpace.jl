# B-spline space

abstract type AbstractFunctionSpace{T} end
@doc raw"""
Abstract type for B-spline space (piecewise polynomial space).

# Examples
```jldoctest
julia> BSplineSpace <: AbstractBSplineSpace
true

julia> UniformBSplineSpace <: AbstractBSplineSpace
true
```
"""
abstract type AbstractBSplineSpace{p,T} <: AbstractFunctionSpace{T} end

# Broadcast like a scalar
Base.Broadcast.broadcastable(P::AbstractBSplineSpace) = Ref(P)

# Equality
@inline Base.:(==)(P1::AbstractBSplineSpace{p}, P2::AbstractBSplineSpace{p}) where p = knotvector(P1) == knotvector(P2)
@inline Base.:(==)(P1::AbstractBSplineSpace{p1}, P2::AbstractBSplineSpace{p2}) where {p1, p2} = false

@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```

# Examples
```jldoctest
julia> p = 2
2

julia> k = KnotVector([1,3,5,6,8,9])
KnotVector([1, 3, 5, 6, 8, 9])

julia> BSplineSpace{p}(k)
BSplineSpace{2, Int64}(KnotVector([1, 3, 5, 6, 8, 9]))
```
"""
struct BSplineSpace{p, T<:Real} <: AbstractBSplineSpace{p,T}
    knotvector::KnotVector{T}
    global unsafe_bsplinespace(::Val{p}, k::AbstractKnotVector{T}) where {p,T} = new{p,T}(k)
end
function BSplineSpace{p}(k::AbstractKnotVector) where p
    if p < 0
        throw(DomainError(p, "degree of polynominal must be non-negative"))
    end
    unsafe_bsplinespace(Val{p}(), k)
end

"""
Convert AbstractBSplineSpace to BSplineSpace
"""
function BSplineSpace(P::AbstractBSplineSpace{p}) where p
    return BSplineSpace{p}(knotvector(P))
end
function BSplineSpace{p}(P::AbstractBSplineSpace{p}) where p
    return BSplineSpace{p}(knotvector(P))
end
function BSplineSpace{p,T}(P::AbstractBSplineSpace{p}) where {p,T}
    return unsafe_bsplinespace(Val{p}(),KnotVector{T}(knotvector(P)))
end

bsplinespace(P::AbstractBSplineSpace) = P

@inline function degree(::AbstractBSplineSpace{p}) where p
    return p
end

@inline function knotvector(P::BSplineSpace)
    return P.knotvector
end

function domain(P::AbstractBSplineSpace)
    p = degree(P)
    k = knotvector(P)
    return k[1+p]..k[end-p]
end

@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\# k - p -1
```
"""
function dim(bsplinespace::AbstractBSplineSpace{p}) where p
    k = knotvector(bsplinespace)
    return length(k) - p - 1
end

@doc raw"""
Check inclusive relationship between B-spline spaces.
```math
\mathcal{P}[p,k]
\subseteq\mathcal{P}[p',k']
```

# Examples
```jldoctest
julia> P1 = BSplineSpace{1}(KnotVector([1,3,5,8]));

julia> P2 = BSplineSpace{1}(KnotVector([1,3,5,6,8,9]));

julia> P3 = BSplineSpace{2}(KnotVector([1,1,3,3,5,5,8,8]));

julia> P1 ⊆ P2
true

julia> P1 ⊆ P3
true

julia> P2 ⊆ P3
false

julia> P2 ⊈ P3
true
```
"""
function Base.issubset(P::AbstractBSplineSpace{p}, P′::AbstractBSplineSpace{p′}) where {p, p′}
    k = knotvector(P)
    k′ = knotvector(P′)
    p₊ = p′ - p

    return p₊ ≥ 0 && (k + p₊ * unique(k) ⊆ k′)
end

@doc raw"""
Check inclusive relationship between B-spline spaces.
```math
\mathcal{P}[p,k]
\sqsubseteq\mathcal{P}[p',k']
\Leftrightarrow
\mathcal{P}[p,k]|_{[k_{p+1},k_{l-p}]}
\subseteq\mathcal{P}[p',k']|_{[k'_{p'+1},k'_{l'-p'}]}
```
"""
function issqsubset(P::AbstractBSplineSpace{p}, P′::AbstractBSplineSpace{p′}) where {p, p′}
    k = knotvector(P)
    k′ = knotvector(P′)
    p₊ = p′ - p

    if p₊ < 0
        return false
    elseif domain(P) ≠ domain(P′)
        return false
    end

    inner_knotvector = k[p+2:end-p-1]
    inner_knotvector′ = k′[p′+2:end-p′-1]

    return inner_knotvector + p₊ * unique(inner_knotvector) ⊆ inner_knotvector′
end

const ⊑ = issqsubset
⊒(l, r) = r ⊑ l
⋢(l, r) = !⊑(l, r)
⋣(l, r) = r ⋢ l

≃(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace) = (P1 ⊑ P2) & (P2 ⊑ P1)

Base.:⊊(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊆ B)
Base.:⊋(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊇ B)
⋤(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊑ B)
⋥(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊒ B)

function isdegenerate_R(P::AbstractBSplineSpace{p}, i::Integer) where p
    k = knotvector(P)
    return k[i] == k[i+p+1]
end

function isdegenerate_I(P::AbstractBSplineSpace{p}, i::Integer) where p
    return iszero(width(bsplinesupport(P,i) ∩ domain(P)))
end

function _iszeros(P::AbstractBSplineSpace{p}) where p
    return [isdegenerate_R(P,i) for i in 1:dim(P)]
end

@doc raw"""
Check if given B-spline space is non-degenerate.

# Examples
```jldoctest
julia> isnondegenerate(BSplineSpace{2}(KnotVector([1,3,5,6,8,9])))
true

julia> isnondegenerate(BSplineSpace{1}(KnotVector([1,3,3,3,8,9])))
false
```
"""
isnondegenerate(P::AbstractBSplineSpace)

@doc raw"""
Check if given B-spline space is degenerate.

# Examples
```jldoctest
julia> isdegenerate(BSplineSpace{2}(KnotVector([1,3,5,6,8,9])))
false

julia> isdegenerate(BSplineSpace{1}(KnotVector([1,3,3,3,8,9])))
true
```
"""
isdegenerate(P::AbstractBSplineSpace)

for (f, fnon) in ((:isdegenerate_R, :isnondegenerate_R), (:isdegenerate_I, :isnondegenerate_I))
    @eval function $f(P::AbstractBSplineSpace)
        for i in 1:dim(P)
            $f(P,i) && return true
        end
        return false
    end
    @eval $fnon(P::AbstractBSplineSpace, i::Int) = !$f(P, i)
    @eval $fnon(P::AbstractBSplineSpace) = !$f(P)
end
isdegenerate(P::AbstractBSplineSpace, i::Integer) = isdegenerate_R(P,i)
isdegenerate(P::AbstractBSplineSpace) = isdegenerate_R(P)
isnondegenerate(P::AbstractBSplineSpace, i::Integer) = isnondegenerate_R(P, i)
isnondegenerate(P::AbstractBSplineSpace) = isnondegenerate_R(P)

"""
Exact dimension of a B-spline space.

# Examples
```jldoctest
julia> exactdim(BSplineSpace{1}(KnotVector([1,2,3,4,5])))
3

julia> exactdim(BSplineSpace{1}(KnotVector([1,2,2,2,4])))
2
```
"""
function exactdim(P::AbstractBSplineSpace)
    n = dim(P)
    for i in 1:dim(P)
        n -= isdegenerate_R(P,i)
    end
    return n
end

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

julia> P = BSplineSpace{2}(k)
BSplineSpace{2, Float64}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> bsplinesupport(P,1)
0.0..5.5

julia> bsplinesupport(P,2)
1.5..8.0
```
"""
bsplinesupport(P::AbstractBSplineSpace, i::Integer) = bsplinesupport_R(P,i)

function bsplinesupport_R(P::AbstractBSplineSpace{p}, i::Integer) where p
    k = knotvector(P)
    return k[i]..k[i+p+1]
end

function bsplinesupport_I(P::AbstractBSplineSpace{p}, i::Integer) where p
    return bsplinesupport_R(P,i) ∩ domain(P)
end

@doc raw"""
Internal methods for obtaining a B-spline space with one degree lower.
```math
\begin{aligned}
\mathcal{P}[p,k] &\mapsto \mathcal{P}[p-1,k] \\
D^r\mathcal{P}[p,k] &\mapsto D^{r-1}\mathcal{P}[p-1,k]
\end{aligned}
```
"""
_lower

_lower(P::BSplineSpace{p,T}) where {p,T} = BSplineSpace{p-1}(knotvector(P))

"""
Return an index of a interval in the domain of B-spline space

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]);

julia> P = BSplineSpace{2}(k);

julia> domain(P)
2.5..9.0

julia> intervalindex(P,2.6)
1

julia> intervalindex(P,5.6)
2

julia> intervalindex(P,8.5)
3

julia> intervalindex(P,9.5)
3
```
"""
function intervalindex(P::AbstractBSplineSpace{p},t::Real) where p
    k = knotvector(P)
    l = length(k)
    v = view(_vec(k),2+p:l-p-1)
    return searchsortedlast(v,t)+1
end

"""
Expand B-spline space with given additional degree and knotvector.
"""
function expandspace_I(P::BSplineSpace{p,T}; p₊::Integer=0, k₊::KnotVector=KnotVector{T}()) where {p,T}
    k = knotvector(P)
    k̂ = unique(k[1+p:end-p])
    p′ = p + p₊
    k′ = k + p₊*k̂ + k₊
    P′ = BSplineSpace{p′}(k′)
    return P′
end

"""
Expand B-spline space with given additional degree and knotvector.
"""
function expandspace_R(P::BSplineSpace{p,T}; p₊::Integer=0, k₊::KnotVector=KnotVector{T}()) where {p,T}
    k = knotvector(P)
    p′ = p + p₊
    k′ = k + p₊*k
    P′ = BSplineSpace{p′}(k′)
    return P′
end

"""
Expand B-spline space with given additional degree and knotvector.
"""
function expandspace(P::BSplineSpace{p,T}; p₊=0, k₊=KnotVector{T}()) where {p,T}
    expandspace_I(P,p₊=p₊,k₊=k₊)
end
