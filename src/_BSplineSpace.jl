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
    return length(k)-p-1
end

"""
Check inclusive relationship between B-spline spaces.
"""
function Base.:⊆(P::AbstractBSplineSpace, P′::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    p′ = degree(P′)
    k′ = knots(P′)
    p₊ = p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function iszeros(P::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    n = dim(P)
    return [k[i] == k[i+p+1] for i ∈ 1:n]
end

function isproper(P::AbstractBSplineSpace)
    return !|(iszeros(P)...)
end

function properdim(P::AbstractBSplineSpace)
    return dim(P) - sum(iszeros(P))
end
