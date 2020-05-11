# B-spline space
struct FastBSplineSpace{T}
    vector::Array{Float64,1}
    function FastBSplineSpace(p::Int, knots::Knots)
        if p < 0
            error("degree of polynominal must be non-negative")
        end
        new{p}(knots.vector)
    end
end

const f𝒫 = FastBSplineSpace

function Knots(P::FastBSplineSpace)
    return Knots(P.vector)
end


@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\sharp k - p -1
```
"""
function dim(bsplinespace::BSplineSpace{p})
    k=bsplinespace.vector
    return length(k)-p-1
end

"""
Check inclusive relationship between B-spline spaces.
"""
function Base.:⊆(P::FastBSplineSpace{p}, P′::FastBSplineSpace{p′})
    k = Knots(P)
    k′ = Knots(P′)
    p₊ = p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function iszeros(P::FastBSplineSpace{p})
    k = P.vector
    n = dim(P)
    return [k[i] == k[i+p+1] for i ∈ 1:n]
end

function isproper(P::FastBSplineSpace)
    return !|(iszeros(P)...)
end

function bsplinesupport(P::FastBSplineSpace{p})
    k = P.vector
    return [k[i]..k[i+p+1] for i ∈ 1:dim(P)]
end

function properdim(P::BSplineSpace)
    return dim(P) - sum(iszeros(P))
end
