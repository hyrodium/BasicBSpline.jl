# B-spline space

abstract type AbstractBSplineSpace end

@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```
"""
struct BSplineSpace <: AbstractBSplineSpace
    degree::Int
    knots::Knots
    function BSplineSpace(degree::Integer, knots::Knots)
        if degree < 0
            error("degree of polynominal must be non-negative")
        end
        new(degree, knots)
    end
end

"""
convert AbstractBSplineSpace to BSplineSpace
"""
function BSplineSpace(P::AbstractBSplineSpace)
    return BSplineSpace(degree(P), knots(P))
end

function degree(P::BSplineSpace)
    return P.degree
end

function knots(P::BSplineSpace)
    return P.knots
end

function bsplineunity(P::AbstractBSplineSpace)
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
function dim(bsplinespace::AbstractBSplineSpace)
    p = degree(bsplinespace)
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
function Base.:⊆(P::AbstractBSplineSpace, P′::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    p′ = degree(P′)
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
function issqsubset(P::AbstractBSplineSpace, P′::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    p′ = degree(P′)
    k′ = knots(P′)
    p₊ = p′ - p

    if p₊ < 0
        return false
    elseif bsplineunity(P) ≠ bsplineunity(P′)
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

≃(P1, P2) = (P1 ⊑ P2)&(P2 ⊑ P1)

function iszeros(P::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    n = dim(P)
    return [k[i] == k[i+p+1] for i in 1:n]
end

function isproper(P::AbstractBSplineSpace)
    return !|(iszeros(P)...)
end

function properdim(P::AbstractBSplineSpace)
    return dim(P) - sum(iszeros(P))
end
