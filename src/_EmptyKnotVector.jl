struct EmptyKnotVector{T} <: AbstractKnotVector{T} end

EmptyKnotVector() = EmptyKnotVector{Bool}()

# + AbstractKnotVector
Base.:+(k::AbstractKnotVector{T}, ::EmptyKnotVector{T}) where T<:Real = k
Base.:+(k::AbstractKnotVector{T1}, ::EmptyKnotVector{T2}) where {T1<:Real, T2<:Real} = AbstractKnotVector{promote_type(T1, T2)}(k)
# + EmptyKnotVector
Base.:+(::EmptyKnotVector{T}, ::EmptyKnotVector{T}) where T<:Real = EmptyKnotVector{T}()
Base.:+(::EmptyKnotVector{T1}, ::EmptyKnotVector{T2}) where {T1<:Real, T2<:Real} = EmptyKnotVector{promote_type(T1, T2)}()
# + swap
Base.:+(k1::EmptyKnotVector, k2::AbstractKnotVector) = k2 + k1

function Base.show(io::IO, ::T) where T<:EmptyKnotVector
    print(io, "$(T)()")
end

Base.zero(::Type{<:EmptyKnotVector{T}}) where {T<:Real} = EmptyKnotVector{T}()
Base.zero(::Type{EmptyKnotVector}) = EmptyKnotVector()
Base.zero(::T) where {T<:EmptyKnotVector} = zero(T)
