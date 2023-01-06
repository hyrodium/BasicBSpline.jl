# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    iszero(b) && return zero(U)
    return U(a/b)
end
@inline _d(a,b) = _d(promote(a,b)...)

@doc raw"""
Calculate ``r``-nomial coefficient

    r_nomial(n, k, r)

```math
(1+x+\cdots+x^r)^n = \sum_{k} a_{n,k,r} x^k
```
"""
function r_nomial(n::T,k::T,r::T) where T<:Integer
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
        return r_nomial(n,(r-one(T))*n-k,r)
    elseif n == 1
        return one(T)
    elseif n == 2
        return k+one(T)
    else
        m = r_nomial(n-one(T),k,r)
        for i in 1:r-1
            m += r_nomial(n-one(T),k-i,r)
        end
        return m
    end
end
