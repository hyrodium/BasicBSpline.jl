# Faster B-spline space

@doc raw"""
B-spline space for lower polynomial degree
```math
\mathcal{P}[p,k]
```
This type `FastBSplineSpace` is faster than `BSplineSpace`, but the degree must be equal or less than `MAX_DEGREE`.
"""
struct FastBSplineSpace{p} <: AbstractBSplineSpace{p}
    knots::Knots
    function FastBSplineSpace(p::Integer, knots::Knots)
        if p < 0
            throw(DomainError(p, "degree of polynominal must be non-negative"))
        elseif p > MAX_DEGREE
            throw(DomainError(p, "FastBSpline supports only degree 0,...,$(MAX_DEGREE)"))
        end
        new{p}(knots)
    end
end
function FastBSplineSpace{q}(p::Integer, knots::Knots) where {q}
    FastBSplineSpace(p,knots)
end

"""
convert AbstractBSplineSpace to FastBSplineSpace
"""
function FastBSplineSpace(P::AbstractBSplineSpace)
    return FastBSplineSpace(degree(P), knots(P))
end

# """
# convert AbstractBSplineSpace to FastBSplineSpace
# """
# function Base.convert(::Type{<:FastBSplineSpace}, P::AbstractBSplineSpace)
#     return FastBSplineSpace(degree(P), knots(P))
# end

function degree(::FastBSplineSpace{p}) where {p}
    return p
end

function knots(P::FastBSplineSpace)
    return P.knots
end
