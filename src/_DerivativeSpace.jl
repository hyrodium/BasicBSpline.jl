# Space of derivative of B-spline basis function

struct BSplineDerivativeSpace{r, T<:AbstractBSplineSpace}
    bsplinespace::T
    function BSplineDerivativeSpace{r}(P::T) where {r, T<:AbstractBSplineSpace}
        new{r,T}(P)
    end
end

bsplinespace(dP::BSplineDerivativeSpace) = dP.bsplinespace
knotvector(dP::BSplineDerivativeSpace) = knotvector(bsplinespace(dP))
degree(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = p - r
dim(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = dim(bsplinespace(dP))
properdim(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = properdim(bsplinespace(dP)) - r
intervalindex(dP::BSplineDerivativeSpace,t::Real) = intervalindex(bsplinespace(dP),t)
domain(dP::BSplineDerivativeSpace) = domain(bsplinespace(dP))
_lower(dP::BSplineDerivativeSpace{r}) where r = BSplineDerivativeSpace{r-1}(_lower(bsplinespace(dP)))

function Base.issubset(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}, P′::AbstractBSplineSpace) where {r,p}
    k = knotvector(dP)
    P = BSplineSpace{p-r}(k)
    return P ⊆ P′
end
function Base.issubset(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}, dP′::BSplineDerivativeSpace{0}) where {r,p}
    P′ = bsplinespace(dP′)
    return dP ⊆ P′
end
function Base.issubset(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}, dP′::BSplineDerivativeSpace{r′,<:AbstractBSplineSpace{p}}) where {r,p,r′,p′}
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
        # This might not be correct.
        return false
    end
end
function Base.issubset(P::AbstractBSplineSpace, dP′::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p}
    # This might not be correct.
    return false
end

# TODO: Add issqsubset
