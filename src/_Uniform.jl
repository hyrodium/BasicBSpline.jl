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

@inline function uniform_bsplinebasisall_kernel(::Val{0},t)
    SVector{1}(one(t))
end
@inline function uniform_bsplinebasisall_kernel(::Val{1},t)
    SVector{2}(1-t,t)
end
@generated function uniform_bsplinebasisall_kernel(::Val{p}, t) where p
    bs = [Symbol(:b,i) for i in 1:p]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    K1s = [:(($(j)-t)/$(p)) for j in 1:p]
    K2s = [:((t-$(j-p))/$(p)) for j in 1:p]
    b = Expr(:tuple, bs...)
    B = Expr(:tuple, Bs...)
    exs = [:($(Bs[j+1]) = ($(K1s[j+1])*$(bs[j+1]) + $(K2s[j])*$(bs[j]))) for j in 1:p-1]
    Expr(
        :block,
        :($(Expr(:meta, :inline))),
        :($b = uniform_bsplinebasisall_kernel(Val{$(p-1)}(),t)),
        :($(Bs[1]) = $(K1s[1])*$(bs[1])),
        exs...,
        :($(Bs[p+1]) = $(K2s[p])*$(bs[p])),
        :(return SVector($(B)))
    )
end
