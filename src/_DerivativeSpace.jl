# Space of derivative of B-spline basis function
@doc raw"""
    BSplineDerivativeSpace{r}(P::BSplineSpace)

Construct derivative of B-spline space from given differential order and B-spline space.
```math
D^{r}(\mathcal{P}[p,k])
=\left\{t \mapsto \left. \frac{d^r f}{dt^r}(t) \  \right| \ f \in \mathcal{P}[p,k] \right\}
```

# Examples
```jldoctest
julia> P = BSplineSpace{2}(KnotVector([1,2,3,4,5,6]))
BSplineSpace{2, Int64}(KnotVector([1, 2, 3, 4, 5, 6]))

julia> dP = BSplineDerivativeSpace{1}(P)
BSplineDerivativeSpace{1, BSplineSpace{2, Int64}, Int64}(BSplineSpace{2, Int64}(KnotVector([1, 2, 3, 4, 5, 6])))

julia> degree(P), degree(dP)
(2, 1)
```
"""
struct BSplineDerivativeSpace{r, S<:BSplineSpace, T} <: AbstractFunctionSpace{T}
    bsplinespace::S
    function BSplineDerivativeSpace{r,S,T}(P::S) where {r, S<:BSplineSpace{p,T}} where {p,T}
        new{r,S,T}(P)
    end
end

function BSplineDerivativeSpace{r}(P::S) where {r, S<:BSplineSpace{p,T}} where {p,T}
    BSplineDerivativeSpace{r,S,T}(P)
end
function BSplineDerivativeSpace{r,S}(P::S) where {r, S<:BSplineSpace{p,T}} where {p,T}
    BSplineDerivativeSpace{r,S,T}(P)
end
function BSplineDerivativeSpace{r,S}(dP::BSplineDerivativeSpace{r,S}) where {r,S}
    dP
end
function BSplineDerivativeSpace{r,S}(dP::BSplineDerivativeSpace{r}) where {r,S}
    BSplineDerivativeSpace{r,S}(S(bsplinespace(dP)))
end

# Broadcast like a scalar
Base.Broadcast.broadcastable(dP::BSplineDerivativeSpace) = Ref(dP)

# Equality
@inline Base.:(==)(dP1::BSplineDerivativeSpace{r}, dP2::BSplineDerivativeSpace{r}) where r = bsplinespace(dP1) == bsplinespace(dP2)
@inline Base.:(==)(dP1::BSplineDerivativeSpace{r1}, dP2::BSplineDerivativeSpace{r2}) where {r1, r2} = false

bsplinespace(dP::BSplineDerivativeSpace) = dP.bsplinespace
knotvector(dP::BSplineDerivativeSpace) = knotvector(bsplinespace(dP))
degree(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}) where {r,p} = p - r
dim(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}) where {r,p} = dim(bsplinespace(dP))
exactdim(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}) where {r,p} = exactdim(bsplinespace(dP)) - r
intervalindex(dP::BSplineDerivativeSpace,t::Real) = intervalindex(bsplinespace(dP),t)
domain(dP::BSplineDerivativeSpace) = domain(bsplinespace(dP))
_lower(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r-1}(_lower(bsplinespace(dP)))

@doc raw"""
    derivative(::BSplineDerivativeSpace{r}) -> BSplineDerivativeSpace{r+1}
    derivative(::BSplineSpace) -> BSplineDerivativeSpace{1}

Derivative of B-spline related space.

# Examples
```jldoctest
julia> BSplineSpace{2}(KnotVector(0:5))
BSplineSpace{2, Int64}(KnotVector([0, 1, 2, 3, 4, 5]))

julia> BasicBSpline.derivative(ans)
BSplineDerivativeSpace{1, BSplineSpace{2, Int64}, Int64}(BSplineSpace{2, Int64}(KnotVector([0, 1, 2, 3, 4, 5])))

julia> BasicBSpline.derivative(ans)
BSplineDerivativeSpace{2, BSplineSpace{2, Int64}, Int64}(BSplineSpace{2, Int64}(KnotVector([0, 1, 2, 3, 4, 5])))
```
"""
derivative(P::BSplineSpace) = BSplineDerivativeSpace{1}(P)
derivative(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r+1}(bsplinespace(dP))

function Base.issubset(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}, P′::BSplineSpace) where {r,p}
    k = knotvector(dP)
    P = BSplineSpace{p-r}(k)
    return P ⊆ P′
end
function Base.issubset(dP::BSplineDerivativeSpace, dP′::BSplineDerivativeSpace{0})
    P′ = bsplinespace(dP′)
    return dP ⊆ P′
end
function Base.issubset(dP::BSplineDerivativeSpace{r}, dP′::BSplineDerivativeSpace{r′}) where {r,r′}
    if r > r′
        P = bsplinespace(dP)
        P′ = bsplinespace(dP′)
        _dP = BSplineDerivativeSpace{r-r′}(P)
        return _dP ⊆ P′
    elseif r == r′
        P = bsplinespace(dP)
        P′ = bsplinespace(dP′)
        return P ⊆ P′
    else
        return false
    end
end
function Base.issubset(P::BSplineSpace, dP′::BSplineDerivativeSpace{0})
    P′ = bsplinespace(dP′)
    return P ⊆ P′
end
function Base.issubset(P::BSplineSpace, dP′::BSplineDerivativeSpace)
    return false
end

# TODO: Add issqsubset
