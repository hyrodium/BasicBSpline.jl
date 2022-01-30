# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    return ifelse(iszero(b), zero(U), U(a/b))
end
@inline _d(a,b) = _d(promote(a,b)...)

# Create ininity from given type
@inline _inf(::Type{T}) where T = one(T)/zero(T)
