# B-spline space

abstract type AbstractBSplineSpace{p,T} end

@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```
This type `BSplineSpace` is slower than `FastBSplineSpace`, but this type is not limited with degree.
"""
struct BSplineSpace{p, T<:Real} <: AbstractBSplineSpace{p,T}
    knots::Knots{T}
    global unsafe_bsplinespace(::Val{p}, k::Knots{T}) where {p,T} = new{p,T}(k)
end
function BSplineSpace{p}(k::Knots) where p
    # TOOD: add error handling like this:
    # throw(DomainError(p, "degree of polynominal must be non-negative"))
    unsafe_bsplinespace(Val{p}(), k)
end

"""
convert AbstractBSplineSpace to BSplineSpace
"""
function BSplineSpace(P::AbstractBSplineSpace{p}) where p
    return BSplineSpace{degree(P)}(knots(P))
end

@inline function degree(::BSplineSpace{p}) where p
    return p
end

@inline function knots(P::BSplineSpace)
    return P.knots
end

function domain(P::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
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
    k = knots(bsplinespace)
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
    k = knots(P)
    k′ = knots(P′)
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
    k = knots(P)
    k′ = knots(P′)
    p₊ = p′ - p

    if p₊ < 0
        return false
    elseif domain(P) ≠ domain(P′)
        return false
    end

    inner_knots = k[p+2:end-p-1]
    inner_knots′ = k′[p′+2:end-p′-1]

    return inner_knots + p₊ * unique(inner_knots) ⊆ inner_knots′
end

const ⊑ = issqsubset
⊒(l, r) = r ⊑ l
⋢(l, r) = !⊑(l, r)
⋣(l, r) = r ⋢ l

≃(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace) = (P1 ⊑ P2) & (P2 ⊑ P1)

function iszeros(P::AbstractBSplineSpace{p}) where p
    k = knots(P)
    n = dim(P)
    return [k[i] == k[i+p+1] for i in 1:n]
end

function isproper(P::AbstractBSplineSpace)
    return !any(iszeros(P))
end

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
    k = knots(P)
    return k[i]..k[i+p+1]
end

function bsplinesupport(P::AbstractBSplineSpace{p}) where p
    k = knots(P)
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
_lower(P::BSplineSpace{p,T}) where {p,T} = BSplineSpace{p-1}(knots(P))

"""
TODO: add docstring
"""
function intervalindex(P::AbstractBSplineSpace{p},t::Real) where p
    k = knots(P)
    l = length(k)
    v = view(k.vector,2+p:l-p-1)
    return searchsortedlast(v,t)+1
end

"""
Expand B-spline space with given additional degree and knots.
"""
function expandspace(P::BSplineSpace{p,T}; p₊::Integer=0, k₊::Knots{T}=Knots{T}()) where {p,T}
    k = knots(P)
    k0 = unique(k[1+p:end-p])
    p′ = p + p₊
    k′ = k + p₊*k0 + k₊
    P′ = BSplineSpace{p′}(k′)
    return P′
end
