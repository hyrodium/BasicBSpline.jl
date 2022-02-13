# Uniform BSpline Basis

## Kernel funcitons of uniform B-spline basis function
@inline function uniform_bsplinebasisall_kernel(::Val{0},t::T) where T<:Real
    U = StaticArrays.arithmetic_closure(T)
    SVector{1,U}(one(t))
end
@inline function uniform_bsplinebasisall_kernel(::Val{1},t::T) where T<:Real
    U = StaticArrays.arithmetic_closure(T)
    SVector{2,U}(1-t,t)
end
@generated function uniform_bsplinebasisall_kernel(::Val{p}, t::T) where {p, T<:Real}
    bs = [Symbol(:b,i) for i in 1:p]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    K1s = [:(($(j)-t)/$(p)) for j in 1:p]
    K2s = [:((t-$(j-p))/$(p)) for j in 1:p]
    b = Expr(:tuple, bs...)
    B = Expr(:tuple, Bs...)
    exs = [:($(Bs[j+1]) = ($(K1s[j+1])*$(bs[j+1]) + $(K2s[j])*$(bs[j]))) for j in 1:p-1]
    Expr(
        :block,
        :($(Expr(:meta, :inline))),
        :($b = uniform_bsplinebasisall_kernel(Val{$(p-1)}(),t)),
        :($(Bs[1]) = $(K1s[1])*$(bs[1])),
        exs...,
        :($(Bs[p+1]) = $(K2s[p])*$(bs[p])),
        :(return SVector($(B)))
    )
end

@inline function uniform_bsplinebasis_kernel(::Val{0},t::T) where T<:Real
    U = StaticArrays.arithmetic_closure(T)
    return U(zero(t) ≤ t < one(t))
end
@inline function uniform_bsplinebasis_kernel(::Val{1},t::T) where T<:Real
    U = StaticArrays.arithmetic_closure(T)
    return U(max(1-abs(t-1), zero(t)))
end
@inline function uniform_bsplinebasis_kernel(::Val{p},t::T) where {p,T<:Real}
    return (t*uniform_bsplinebasis_kernel(Val{p-1}(),t) + (p+1-t)*uniform_bsplinebasis_kernel(Val{p-1}(),t-1))/p
end

## bsplinebasis
@inline function bsplinebasis(P::UniformBSplineSpace{p,T,R},i::Integer,t::S) where {p, T, R<:AbstractUnitRange, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (1 ≤ i ≤ dim(P)) || throw(DomainError(i, "index of B-spline basis function is out of range."))
    return U(uniform_bsplinebasis_kernel(Val{p}(),t-i+1-k[1]))
end
@inline function bsplinebasis(P::UniformBSplineSpace{p,T},i::Integer,t::S) where {p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (1 ≤ i ≤ dim(P)) || throw(DomainError(i, "index of B-spline basis function is out of range."))
    a = @inbounds k.vector[i]
    b = @inbounds k.vector[i+p+1]
    return U(uniform_bsplinebasis_kernel(Val{p}(),(p+1)*(t-a)/(b-a)))
end
@inline bsplinebasis₋₀(P::UniformBSplineSpace,i::Integer,t::Real) = bsplinebasis(P,i,t)
@inline bsplinebasis₊₀(P::UniformBSplineSpace,i::Integer,t::Real) = bsplinebasis(P,i,t)

## bsplinebasisall
@inline function bsplinebasisall(P::UniformBSplineSpace{p,T,R},i::Integer,t::S) where {p, T, R<:AbstractUnitRange, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (0 ≤ i ≤ length(k)-2p) || throw(DomainError(i, "index of interval is out of range."))
    return uniform_bsplinebasisall_kernel(Val{p}(),U(t-i+1-p-k[1]))
end
@inline function bsplinebasisall(P::UniformBSplineSpace{p,T},i::Integer,t::S) where {p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (0 ≤ i ≤ length(k)-2p) || throw(DomainError(i, "index of interval is out of range."))
    a = @inbounds k.vector[i+p]
    b = @inbounds k.vector[i+p+1]
    return uniform_bsplinebasisall_kernel(Val{p}(),U((t-a)/(b-a)))
end
