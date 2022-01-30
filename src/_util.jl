# Division, devided by zero will be zero
@inline function _d(a::T,b::T) where T
    U=StaticArrays.arithmetic_closure(T)
    iszero(b) && return zero(U)
    return U(a/b)
end
@inline _d(a,b) = _d(promote(a,b)...)
