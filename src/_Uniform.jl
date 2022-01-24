struct UniformKnotVector{T,R<:StepRangeLen} <: AbstractKnotVector{T}
    vector::R
    global unsafe_uniformknotvector(v::R) where R<:StepRangeLen{T} where T = new{T,R}(v)
end
UniformKnotVector(v::AbstractRange) = unsafe_uniformknotvector(StepRangeLen(sort(v)))

Base.in(r::Real, k::UniformKnotVector) = in(r, k.vector)
Base.getindex(k::UniformKnotVector, i::Integer) = k.vector[i]
Base.getindex(k::UniformKnotVector, v::AbstractVector{<:Integer}) = KnotVector(k.vector[v])
Base.getindex(k::UniformKnotVector, v::AbstractRange{<:Integer}) = UniformKnotVector(sort(k.vector[v]))

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
