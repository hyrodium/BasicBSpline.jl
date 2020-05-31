# B-spline space
struct FastBSplineSpace{T} <: AbstractBSplineSpace
    vector::Array{Float64,1}
    function FastBSplineSpace(p::Int, knots::Knots)
        if p < 0
            error("degree of polynominal must be non-negative")
        elseif p > MAX_DEGREE
            error("FastBSpline supports only degree 0 , ... , 3")
        end
        new{p}(knots.vector)
    end
end

# function knots(P::FastBSplineSpace)
#     return Knots(P.vector)
# end


@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\sharp k - p -1
```
"""
function dim(P::FastBSplineSpace{p}) where p
    k=P.vector
    return length(k)-p-1
end

@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\sharp k - p -1
```
"""
function BSplineSpace(P::FastBSplineSpace{p}) where p
    BSplineSpace(p, Knots(P.vector))
end


"""
Check inclusive relationship between B-spline spaces.
"""
function Base.:⊆(P::FastBSplineSpace{p}, P′::FastBSplineSpace{p′}) where p where p′
    k = knots(P)
    k′ = knots(P′)
    p₊ = p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function iszeros(P::FastBSplineSpace{p}) where p
    k = P.vector
    n = dim(P)
    return [k[i] == k[i+p+1] for i ∈ 1:n]
end

function isproper(P::FastBSplineSpace)
    return !|(iszeros(P)...)
end

function properdim(P::FastBSplineSpace)
    return dim(P) - sum(iszeros(P))
end

function degree(P::FastBSplineSpace{p}) where p
    return p
end
function knots(P::FastBSplineSpace)
    return Knots(P.vector)
end
