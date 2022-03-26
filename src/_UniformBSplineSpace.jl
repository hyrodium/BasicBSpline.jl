# Uniform BSpline Space

struct UniformBSplineSpace{p, T<:Real, R} <: AbstractBSplineSpace{p,T}
    knotvector::UniformKnotVector{T,R}
    global unsafe_uniformbsplinespace(::Val{p}, k::UniformKnotVector{T,R}) where {p,T,R} = new{p,T,R}(k)
end
function UniformBSplineSpace{p}(k::UniformKnotVector) where p
    if p < 0
        throw(DomainError(p, "degree of polynominal must be non-negative"))
    end
    unsafe_uniformbsplinespace(Val{p}(), k)
end
# TODO: add constructor like this
# function UniformBSplineSpace{p,T,R}(k::K) where {p,K<:UniformKnotVector{T,R}} where {T,R}
#     UniformBSplineSpace{p}()
# end
function UniformBSplineSpace{p,T,R}(P::UniformBSplineSpace{p,T,R}) where {p,R<:AbstractRange{T}} where T
    P
end
function UniformBSplineSpace{p,T,R}(P::AbstractBSplineSpace{p}) where {p,R<:AbstractRange{T}} where T
    k = knotvector(P)
    unsafe_uniformbsplinespace(Val{p}(), UniformKnotVector{T,R}(k))
end

knotvector(P::UniformBSplineSpace) = P.knotvector

_lower(P::UniformBSplineSpace{p,T}) where {p,T} = UniformBSplineSpace{p-1}(knotvector(P))
