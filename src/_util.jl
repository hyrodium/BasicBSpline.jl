# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    iszero(b) && return zero(U)
    return U(a/b)
end
@inline _d(a,b) = _d(promote(a,b)...)

# Euler's triangle number (http://oeis.org/wiki/Eulerian_numbers,_triangle_of)
eulertriangle(n,k) = sum((-1)^j*binomial(n+1,j)*(k-j+1)^n for j in 0:k)

# Left division function
_leftdivision(A,b) = inv(A)*b
_leftdivision(A,b::AbstractVector{<:Number}) = A\b
# TODO: add more methods for left division (e.g. b::Vector{<:SVector})

trinomial(n,k) = sum((-1)^j*binomial(n,j)*binomial(2n-2j,n-k-j) for j in 0:n)
