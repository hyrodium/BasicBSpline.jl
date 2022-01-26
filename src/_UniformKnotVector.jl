struct UniformKnotVector{T,R<:AbstractRange} <: AbstractKnotVector{T}
    vector::R
    global unsafe_uniformknotvector(::Type{T}, v::R) where R<:AbstractRange{T} where T = new{T,R}(v)
end
UniformKnotVector(v::AbstractRange{T}) where T = unsafe_uniformknotvector(T,sort(v))
UniformKnotVector(k::UniformKnotVector) = k
UniformKnotVector{T}(k::UniformKnotVector) where T = unsafe_uniformknotvector(T, k.vector)
UniformKnotVector{T,R}(k::UniformKnotVector) where R<:AbstractRange{T} where T = unsafe_uniformknotvector(T, R(k.vector))

_vec(k::UniformKnotVector) = k.vector

Base.getindex(k::UniformKnotVector, v::AbstractVector{<:Integer}) = KnotVector(k.vector[v])
Base.getindex(k::UniformKnotVector, v::AbstractRange{<:Integer}) = UniformKnotVector(sort(k.vector[v]))

Base.length(k::UniformKnotVector) = length(k.vector)
Base.firstindex(k::UniformKnotVector) = 1
Base.lastindex(k::UniformKnotVector) = length(k)

function Base.convert(::Type{UniformKnotVector{T,R}},k::UniformKnotVector) where {T,R}
    UniformKnotVector{T,R}(k)
end

function Base.promote_rule(::Type{KnotVector{T}}, ::Type{UniformKnotVector{S,R}}) where {T,S,R}
    KnotVector{promote_type(T,S)}
end

Base.unique(k::UniformKnotVector) = UniformKnotVector(unique(k.vector))

Base.collect(k::UniformKnotVector) = collect(k.vector)

Base.:+(k1::UniformKnotVector{T1},k2::UniformKnotVector{T2}) where {T1,T2} = KnotVector{promote_type(T1,T2)}([k1.vector;k2.vector])

Base.zero(::UniformKnotVector{T}) where T = KnotVector{T}()
