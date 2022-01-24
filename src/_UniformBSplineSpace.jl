# Uniform BSpline Space

struct UniformBSplineSpace{p, T<:Real, R} <: AbstractBSplineSpace{p,T}
    knotvector::UniformKnotVector{T,R}
    global unsafe_uniformbsplinespace(::Val{p}, k::UniformKnotVector{T,N}) where {p,T,N} = new{p,T,N}(k)
end
function UniformBSplineSpace{p}(k::UniformKnotVector) where p
    if p < 0
        throw(DomainError(p, "degree of polynominal must be non-negative"))
    end
    unsafe_uniformbsplinespace(Val{p}(), k)
end

BasicBSpline.knotvector(P::UniformBSplineSpace) = P.knotvector

_lower(P::UniformBSplineSpace{p,T}) where {p,T} = UniformBSplineSpace{p-1}(knotvector(P))
