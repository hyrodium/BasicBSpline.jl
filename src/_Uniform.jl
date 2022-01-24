struct UniformKnotVector{T,R<:AbstractRange} <: AbstractKnotVector{T}
    vector::R
    global unsafe_uniformknotvector(v::R) where R<:AbstractRange{T} where T = new{T,R}(v)
end
UniformKnotVector(v::AbstractRange) = unsafe_uniformknotvector(sort(v))
UniformKnotVector(k::UniformKnotVector) = k

KnotVector(k::UniformKnotVector{T}) where T = unsafe_knotvector(T,k.vector)
KnotVector{T}(k::UniformKnotVector{S}) where {T,S} = unsafe_knotvector(T,k.vector)

Base.in(r::Real, k::UniformKnotVector) = in(r, k.vector)
Base.getindex(k::UniformKnotVector, i::Integer) = k.vector[i]
Base.getindex(k::UniformKnotVector, v::AbstractVector{<:Integer}) = KnotVector(k.vector[v])
Base.getindex(k::UniformKnotVector, v::AbstractRange{<:Integer}) = UniformKnotVector(sort(k.vector[v]))

Base.length(k::UniformKnotVector) = length(k.vector)
Base.firstindex(k::UniformKnotVector) = 1
Base.lastindex(k::UniformKnotVector) = length(k)

function Base.convert(::Type{KnotVector},k::UniformKnotVector{T}) where T
    unsafe_knotvector(T,k.vector)
end

function Base.convert(::Type{KnotVector{T}},k::UniformKnotVector) where T
    unsafe_knotvector(T,k.vector)
end

function Base.promote_rule(::Type{KnotVector{T}}, ::Type{UniformKnotVector{S,R}}) where {T,S,R}
    KnotVector{promote_type(T,S)}
end

Base.unique(k::UniformKnotVector) = UniformKnotVector(unique(k.vector))
Base.iterate(k::UniformKnotVector) = iterate(k.vector)
Base.iterate(k::UniformKnotVector, i::Integer) = iterate(k.vector, i)
Base.searchsortedfirst(k::UniformKnotVector,t) = searchsortedfirst(k.vector,t)
Base.searchsortedlast(k::UniformKnotVector,t) = searchsortedlast(k.vector,t)
Base.searchsorted(k::UniformKnotVector,t) = searchsorted(k.vector,t)

Base.collect(k::UniformKnotVector) = collect(k.vector)

Base.:+(k1::UniformKnotVector{T1},k2::UniformKnotVector{T2}) where {T1,T2} = KnotVector{promote_type(T1,T2)}([k1.vector;k2.vector])

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

@inline function uniform_bsplinebasis_kernel(::Val{0},t::T) where T<:Real
    return zero(t) ≤ t < one(t)
end
@inline function uniform_bsplinebasis_kernel(::Val{1},t::T) where T<:Real
    return max(1-abs(t-1), zero(t))
end
@inline function uniform_bsplinebasis_kernel(::Val{p},t::T) where {p,T<:Real}
    return (t*uniform_bsplinebasis_kernel(Val{p-1}(),t) + (p+1-t)*uniform_bsplinebasis_kernel(Val{p-1}(),t-1))/p
end
