# B-spline space

@doc raw"""
B-spline space for lower polynomial degree
"""
struct FastBSplineSpace{T} <: AbstractBSplineSpace
    vector::Array{Float64,1}
    function FastBSplineSpace(p::Int, knots::Knots)
        if p < 0
            error("degree of polynominal must be non-negative")
        elseif p > MAX_DEGREE
            error("FastBSpline supports only degree 0 , ... , $(MAX_DEGREE)")
        end
        new{p}(knots.vector)
    end
end

"""
convert AbstractBSplineSpace to FastBSplineSpace
"""
function FastBSplineSpace(P::AbstractBSplineSpace)
    return FastBSplineSpace(degree(P), knots(P))
end

function degree(P::FastBSplineSpace{p}) where p
    return p
end

function knots(P::FastBSplineSpace)
    return Knots(P.vector)
end

@doc raw"""
Retrun FastBSplineSpace if ≤ MAX_DEGREE, or BSplineSpace if not.
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
