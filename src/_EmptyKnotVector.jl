@doc raw"""
Knot vector with zero-element.
```math
k=()
```
This struct is intended for internal use.

# Examples
```jldoctest
julia> import BasicBSpline.EmptyKnotVector

julia> EmptyKnotVector()
BasicBSpline.EmptyKnotVector{Bool}()

julia> EmptyKnotVector{Float64}()
BasicBSpline.EmptyKnotVector{Float64}()
```
"""
struct EmptyKnotVector{T} <: AbstractKnotVector{T} end

EmptyKnotVector() = EmptyKnotVector{Bool}()

# ==
Base.:(==)(k::AbstractKnotVector, ::EmptyKnotVector) = isempty(k)
Base.:(==)(k::EmptyKnotVector, ::EmptyKnotVector) = true
Base.:(==)(k1::EmptyKnotVector, k2::AbstractKnotVector) = (k2 == k1)

Base.isempty(::EmptyKnotVector) = true

_vec(::EmptyKnotVector{T}) where T = T[]

Base.copy(k::EmptyKnotVector{T}) where T = k

# + AbstractKnotVector
Base.:+(k::AbstractKnotVector{T}, ::EmptyKnotVector{T}) where T<:Real = k
function Base.:+(k::AbstractKnotVector{T1}, ::EmptyKnotVector{T2}) where {T1<:Real, T2<:Real}
    T = promote_type(T1, T2)
    if T == T1
        return k
    else
        return AbstractKnotVector{T}(k)
    end
end
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
Base.zero(::Type{<:AbstractKnotVector{T}}) where {T<:Real} = EmptyKnotVector{T}()
Base.zero(::Type{AbstractKnotVector}) = EmptyKnotVector()

Base.unique(k::EmptyKnotVector) = k
