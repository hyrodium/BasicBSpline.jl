# B-spline space

abstract type AbstractBSplineSpace{p,T} end

# Broadcast like a scalar
Base.Broadcast.broadcastable(P::AbstractBSplineSpace) = Ref(P)

@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```
"""
struct BSplineSpace{p, T<:Real} <: AbstractBSplineSpace{p,T}
    knotvector::KnotVector{T}
    global unsafe_bsplinespace(::Val{p}, k::KnotVector{T}) where {p,T} = new{p,T}(k)
end
function BSplineSpace{p}(k::KnotVector) where p
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

@inline function degree(::BSplineSpace{p}) where p
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
=\sharp k - p -1
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
\subseteq\mathcal{P}[p',k']|_{[\sharp k'_{p'+1},k'_{\sharp k'-p'}]}
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

function iszeros(P::AbstractBSplineSpace{p}) where p
    k = knotvector(P)
    n = dim(P)
    return [k[i] == k[i+p+1] for i in 1:n]
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
function isnondegenerate(P::AbstractBSplineSpace)
    return !isdegenerate(P)
end

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
function isdegenerate(P::AbstractBSplineSpace)
    return any(iszeros(P))
end

# This binding will be removed when releasing v0.5.0
Base.@deprecate_binding isproper isnondegenerate true

function properdim(P::AbstractBSplineSpace)
    return dim(P) - sum(iszeros(P))
end

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```
"""
function bsplinesupport(P::AbstractBSplineSpace{p}, i::Integer) where p
    k = knotvector(P)
    return k[i]..k[i+p+1]
end

function bsplinesupport(P::AbstractBSplineSpace{p}) where p
    k = knotvector(P)
    return [k[i]..k[i+p+1] for i in 1:dim(P)]
end

@doc raw"""
Return a B-spline space of one degree lower.
```math
\mathcal{P}[p,k] \mapsto \mathcal{P}[p-1,k]
```
"""
_lower

# TODO: Consider we really need these methods.
# _lower(::Type{AbstractBSplineSpace{p}}) where p = AbstractBSplineSpace{p-1}
# _lower(::Type{AbstractBSplineSpace{p,T}}) where {p,T} = AbstractBSplineSpace{p-1,T}
# _lower(::Type{BSplineSpace{p}}) where p = BSplineSpace{p-1}
# _lower(::Type{BSplineSpace{p,T}}) where {p,T} = BSplineSpace{p-1,T}
_lower(P::BSplineSpace{p,T}) where {p,T} = BSplineSpace{p-1}(knotvector(P))

"""
TODO: add docstring
"""
function intervalindex(P::AbstractBSplineSpace{p},t::Real) where p
    k = knotvector(P)
    l = length(k)
    v = view(k.vector,2+p:l-p-1)
    return searchsortedlast(v,t)+1
end

"""
Expand B-spline space with given additional degree and knotvector.
"""
function expandspace(P::BSplineSpace{p,T}; p₊::Integer=0, k₊::KnotVector{T}=KnotVector{T}()) where {p,T}
    k = knotvector(P)
    k0 = unique(k[1+p:end-p])
    p′ = p + p₊
    k′ = k + p₊*k0 + k₊
    P′ = BSplineSpace{p′}(k′)
    return P′
end
