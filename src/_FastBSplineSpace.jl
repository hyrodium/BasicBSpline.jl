# Faster B-spline space

@doc raw"""
B-spline space for lower polynomial degree
"""
struct FastBSplineSpace{p} <: AbstractBSplineSpace
    knots::Knots
    function FastBSplineSpace(p::Integer, knots::Knots)
        if p < 0
            error("degree of polynominal must be non-negative")
        elseif p > MAX_DEGREE
            error("FastBSpline supports only degree 0 , ... , $(MAX_DEGREE)")
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

function degree(P::FastBSplineSpace{p}) where {p}
    return p
end

function knots(P::FastBSplineSpace)
    return P.knots
end
