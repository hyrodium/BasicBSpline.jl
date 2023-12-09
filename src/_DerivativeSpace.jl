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
BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([1, 2, 3, 4, 5, 6]))

julia> dP = BSplineDerivativeSpace{1}(P)
BSplineDerivativeSpace{1, BSplineSpace{2, Int64, KnotVector{Int64}}, Int64}(BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([1, 2, 3, 4, 5, 6])))

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
function BSplineDerivativeSpace{r,S1}(P::S2) where {r, S1<:BSplineSpace{p,T1}, S2<:BSplineSpace{p,T2}} where {p,T1,T2}
    BSplineDerivativeSpace{r}(S1(P))
end
function BSplineDerivativeSpace{r,S}(dP::BSplineDerivativeSpace{r,S}) where {r,S}
    dP
end
function BSplineDerivativeSpace{r,S1}(dP::BSplineDerivativeSpace{r,S2}) where {r,S1,S2}
    BSplineDerivativeSpace{r,S1}(S1(bsplinespace(dP)))
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
exactdim_R(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}) where {r,p} = exactdim_R(bsplinespace(dP)) - r
exactdim_I(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}) where {r,p} = exactdim_I(bsplinespace(dP)) - r
intervalindex(dP::BSplineDerivativeSpace,t::Real) = intervalindex(bsplinespace(dP),t)
domain(dP::BSplineDerivativeSpace) = domain(bsplinespace(dP))
_lower_R(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r-1}(_lower_R(bsplinespace(dP)))
_lower_I(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r-1}(_lower_I(bsplinespace(dP)))
_iszeros_R(P::BSplineDerivativeSpace) = _iszeros_R(bsplinespace(P))
_iszeros_I(P::BSplineDerivativeSpace) = _iszeros_I(bsplinespace(P))

function isdegenerate_R(dP::BSplineDerivativeSpace, i::Integer)
    return isdegenerate_R(bsplinespace(dP), i)
end

function isdegenerate_I(dP::BSplineDerivativeSpace, i::Integer)
    return isdegenerate_I(bsplinespace(dP), i)
end

@doc raw"""
    derivative(::BSplineDerivativeSpace{r}) -> BSplineDerivativeSpace{r+1}
    derivative(::BSplineSpace) -> BSplineDerivativeSpace{1}

Derivative of B-spline related space.

# Examples
```jldoctest
julia> BSplineSpace{2}(KnotVector(0:5))
BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([0, 1, 2, 3, 4, 5]))

julia> BasicBSpline.derivative(ans)
BSplineDerivativeSpace{1, BSplineSpace{2, Int64, KnotVector{Int64}}, Int64}(BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([0, 1, 2, 3, 4, 5])))

julia> BasicBSpline.derivative(ans)
BSplineDerivativeSpace{2, BSplineSpace{2, Int64, KnotVector{Int64}}, Int64}(BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([0, 1, 2, 3, 4, 5])))
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

function Base.hash(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p}}, h::UInt) where {r,p}
    P = bsplinespace(dP)
    k = knotvector(P)
    hash(BSplineDerivativeSpace{r,<:BSplineSpace{p}}, hash(_vec(k), h))
end
