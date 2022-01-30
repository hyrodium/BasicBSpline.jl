# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    return ifelse(iszero(b), zero(U), U(a/b))
end
@inline _d(a,b) = _d(promote(a,b)...)

# Create ininity from given type
@inline _inf(::Type{T}) where T<:Number = one(T)/zero(T)
@inline _inf(::T) where T<:Number = _inf(T)

# Euler's triangle number (http://oeis.org/wiki/Eulerian_numbers,_triangle_of)
eulertriangle(n,k) = sum((-1)^j*binomial(n+1,j)*(k-j+1)^n for j in 0:k)
