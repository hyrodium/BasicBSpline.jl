struct UniformKnotVector{T,R<:AbstractRange} <: AbstractKnotVector{T}
    vector::R
    global unsafe_uniformknotvector(v::R) where R<:AbstractRange{T} where T = new{T,R}(v)
end
UniformKnotVector(v::AbstractRange) = unsafe_uniformknotvector(sort(v))
UniformKnotVector(k::UniformKnotVector) = k

_vec(k::UniformKnotVector) = k.vector

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

# TODO: remove collect
Base.collect(k::UniformKnotVector) = collect(k.vector)

Base.:+(k1::UniformKnotVector{T1},k2::UniformKnotVector{T2}) where {T1,T2} = KnotVector{promote_type(T1,T2)}([k1.vector;k2.vector])
