# Uniform BSpline Basis

@inline function uniform_bsplinebasisall_kernel(::Val{0},t)
    SVector{1}(one(t))
end
@inline function uniform_bsplinebasisall_kernel(::Val{1},t)
    SVector{2}(1-t,t)
end
@generated function uniform_bsplinebasisall_kernel(::Val{p}, t) where p
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
    return zero(t) ≤ t < one(t)
end
@inline function uniform_bsplinebasis_kernel(::Val{1},t::T) where T<:Real
    return max(1-abs(t-1), zero(t))
end
@inline function uniform_bsplinebasis_kernel(::Val{p},t::T) where {p,T<:Real}
    return (t*uniform_bsplinebasis_kernel(Val{p-1}(),t) + (p+1-t)*uniform_bsplinebasis_kernel(Val{p-1}(),t-1))/p
end
