# B-spline space
@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```
"""
struct BSplineSpace
    degree::Int
    knots::Knots
    function BSplineSpace(degree::Int, knots::Knots)
        if degree < 0
            error("degree of polynominal must be non-negative")
        end
        new(degree, knots)
    end
end

@doc raw"""
Same as BSplineSpace.
```math
\mathcal{P}[p,k]
```
"""
function 𝒫(p::Int,k::Knots)
    if p ≤ MAX_DEGREE
        FastBSplineSpace(p,k)
    else
        BSplineSpace(p,k)
    end
end

@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\sharp k - p -1
```
"""
function dim(bsplinespace::BSplineSpace)
    p=bsplinespace.degree
    k=bsplinespace.knots
    return ♯(k)-p-1
end

"""
Check inclusive relationship between B-spline spaces.
"""
function Base.:⊆(P::BSplineSpace, P′::BSplineSpace)
    p = P.degree
    k = P.knots
    p′ = P′.degree
    k′ = P′.knots
    p₊ = p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function iszeros(P::BSplineSpace)
    p = P.degree
    k = P.knots
    n = dim(P)
    return [k[i] == k[i+p+1] for i ∈ 1:n]
end

function isproper(P::BSplineSpace)
    return !|(iszeros(P)...)
end

function properdim(P::BSplineSpace)
    return dim(P) - sum(iszeros(P))
end
