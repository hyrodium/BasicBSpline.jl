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

function _multinomial(n::T,k::T,r::T) where T<:Integer
    # n must be non-negative
    # r must be larger or equal to one
    if r == 1
        return T(iszero(k))
    elseif r == 2
        return binomial(n,k)
    elseif k < 0
        return zero(T)
    elseif k == 0
        return one(T)
    elseif k == 1
        return n
    elseif k == 2
        return n*(n+1)÷2
    elseif (r-1)*n < 2k
        return _multinomial(n,(r-one(T))*n-k,r)
    elseif (r-1)*n < k
        return zero(T)
    elseif n == 1
        return one(T)
    elseif n == 2
        return k+one(T)
    else
        m = _multinomial(n-one(T),k,r)
        for i in 1:r-1
            m += _multinomial(n-one(T),k-i,r)
        end
        return m
    end
end
