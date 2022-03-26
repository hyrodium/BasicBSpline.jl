# Space of derivative of B-spline basis function

struct BSplineDerivativeSpace{r, S<:AbstractBSplineSpace, T} <: AbstractFunctionSpace{T}
    bsplinespace::S
    function BSplineDerivativeSpace{r,S,T}(P::S) where {r, S<:AbstractBSplineSpace{p,T}} where {p,T}
        new{r,S,T}(P)
    end
end

function BSplineDerivativeSpace{r}(P::S) where {r, S<:AbstractBSplineSpace{p,T}} where {p,T}
    BSplineDerivativeSpace{r,S,T}(P)
end
function BSplineDerivativeSpace{r,S}(P::S) where {r, S<:AbstractBSplineSpace{p,T}} where {p,T}
    BSplineDerivativeSpace{r,S,T}(P)
end
function BSplineDerivativeSpace{r,S}(dP::BSplineDerivativeSpace{r,S}) where {r, S<:AbstractBSplineSpace{p,T}} where {p,T}
    dP
end
function BSplineDerivativeSpace{r,S}(dP::BSplineDerivativeSpace{r}) where {r,S}
    BSplineDerivativeSpace{r,S}(S(bsplinespace(dP)))
end
function Base.convert(::Type{BSplineDerivativeSpace{r,S}}, dP) where {r,S}
    BSplineDerivativeSpace{r,S}(dP)
end

# Broadcast like a scalar
Base.Broadcast.broadcastable(dP::BSplineDerivativeSpace) = Ref(dP)

# Equality
@inline Base.:(==)(dP1::BSplineDerivativeSpace{r}, dP2::BSplineDerivativeSpace{r}) where r = bsplinespace(dP1) == bsplinespace(dP2)
@inline Base.:(==)(dP1::BSplineDerivativeSpace{r1}, dP2::BSplineDerivativeSpace{r2}) where {r1, r2} = false

bsplinespace(dP::BSplineDerivativeSpace) = dP.bsplinespace
knotvector(dP::BSplineDerivativeSpace) = knotvector(bsplinespace(dP))
degree(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = p - r
dim(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = dim(bsplinespace(dP))
exactdim(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = exactdim(bsplinespace(dP)) - r
intervalindex(dP::BSplineDerivativeSpace,t::Real) = intervalindex(bsplinespace(dP),t)
domain(dP::BSplineDerivativeSpace) = domain(bsplinespace(dP))
_lower(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r-1}(_lower(bsplinespace(dP)))

function Base.issubset(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}, P′::AbstractBSplineSpace) where {r,p}
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
function Base.issubset(P::AbstractBSplineSpace, dP′::BSplineDerivativeSpace{0})
    P′ = bsplinespace(dP′)
    return P ⊆ P′
end
function Base.issubset(P::AbstractBSplineSpace, dP′::BSplineDerivativeSpace)
    return false
end

# TODO: Add issqsubset
