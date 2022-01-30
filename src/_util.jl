# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    iszero(b) && return zero(U)
    return U(a/b)
end
@inline _d(a,b) = _d(promote(a,b)...)

# Euler's triangle number (http://oeis.org/wiki/Eulerian_numbers,_triangle_of)
eulertriangle(n,k) = sum((-1)^j*binomial(n+1,j)*(k-j+1)^n for j in 0:k)
