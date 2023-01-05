# Derivative of B-spline basis function

@generated function bsplinebasis₊₀(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p,T}}, i::Integer, t::S) where {r, p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:($U($(ks[i])≤t<$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one($U),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    C_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])-$(Ks[i+1])*$(Bs[i+1])) for i in 1:n]...)
    if r ≤ p
        exs = Expr[]
        for i in 1:p-r
            push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
        end
        for i in p-r+1:p
            push!(exs, :($(K_l(p+2-i)) = $(L_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(C_r(p+1-i))))
        end
        Expr(:block,
            :(v = knotvector(dP).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero($U))
    end
end

@generated function bsplinebasis₋₀(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p,T}}, i::Integer, t::S) where {r, p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:($U($(ks[i])<t≤$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one($U),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    C_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])-$(Ks[i+1])*$(Bs[i+1])) for i in 1:n]...)
    if r ≤ p
        exs = Expr[]
        for i in 1:p-r
            push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
        end
        for i in p-r+1:p
            push!(exs, :($(K_l(p+2-i)) = $(L_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(C_r(p+1-i))))
        end
        Expr(:block,
            :(v = knotvector(dP).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero($U))
    end
end

@generated function bsplinebasis(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p,T}}, i::Integer, t::S) where {r, p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:($U(($(ks[i])≤t<$(ks[i+1])) || ($(ks[i])<t==$(ks[i+1])==v[end]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one($U),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    C_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])-$(Ks[i+1])*$(Bs[i+1])) for i in 1:n]...)
    if r ≤ p
        exs = Expr[]
        for i in 1:p-r
            push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
        end
        for i in p-r+1:p
            push!(exs, :($(K_l(p+2-i)) = $(L_r(i,p+2-i))))
            push!(exs, :($(B_l(p+1-i)) = $(C_r(p+1-i))))
        end
        Expr(:block,
            :(v = knotvector(dP).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero($U))
    end
end

@doc raw"""
    bsplinebasis′₊₀(::AbstractFunctionSpace, ::Integer, ::Real) -> Real

1st derivative of B-spline basis function.
Right-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
`bsplinebasis′₊₀(P, i, t)` is equivalent to `bsplinebasis₊₀(derivative(P), i, t)`.
"""
bsplinebasis′₊₀


@doc raw"""
    bsplinebasis′₋₀(::AbstractFunctionSpace, ::Integer, ::Real) -> Real

1st derivative of B-spline basis function.
Left-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
`bsplinebasis′₋₀(P, i, t)` is equivalent to `bsplinebasis₋₀(derivative(P), i, t)`.
"""
bsplinebasis′₋₀

@doc raw"""
    bsplinebasis′(::AbstractFunctionSpace, ::Integer, ::Real) -> Real

1st derivative of B-spline basis function.
Modified version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
`bsplinebasis′(P, i, t)` is equivalent to `bsplinebasis(derivative(P), i, t)`.
"""
bsplinebasis′

for suffix in ("", "₋₀", "₊₀")
    fname = Symbol(:bsplinebasis, suffix)
    fname′ = Symbol(:bsplinebasis, "′" ,suffix)
    @eval function $(fname′)(P::AbstractFunctionSpace, i::Integer, t::Real)
        return $(fname)(derivative(P), i, t)
    end
end

@generated function bsplinebasisall(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p,T}}, i::Integer, t::S) where {r, p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    bs = [Symbol(:b,i) for i in 1:p]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    K1s = [:($(p)/(k[i+$(j)]-k[i+$(p+j)])) for j in 1:p]
    K2s = [:($(p)/(k[i+$(p+j)]-k[i+$(j)])) for j in 1:p]
    b = Expr(:tuple, bs...)
    B = Expr(:tuple, Bs...)
    exs = [:($(Bs[j+1]) = ($(K1s[j+1])*$(bs[j+1]) + $(K2s[j])*$(bs[j]))) for j in 1:p-1]
    if r ≤ p
        Expr(:block,
            :($(Expr(:meta, :inline))),
            :(k = knotvector(dP)),
            :($b = bsplinebasisall(_lower(dP),i+1,t)),
            :($(Bs[1]) = $(K1s[1])*$(bs[1])),
            exs...,
            :($(Bs[p+1]) = $(K2s[p])*$(bs[p])),
            :(return SVector($(B)))
        )
    else
        Z = Expr(:tuple, [:(zero($U)) for i in 1:p+1]...)
        :(return SVector($(Z)))
    end
end

@inline function bsplinebasisall(dP::BSplineDerivativeSpace{0,<:BSplineSpace{p,T}}, i::Integer, t::Real) where {p, T}
    P = bsplinespace(dP)
    bsplinebasisall(P,i,t)
end

# TODO: add methods for UniformBSplineSpace (and OpenUniformBSplineSpace)
# function bsplinebasis(dP::BSplineDerivativeSpace{r,UniformBSplineSpace{p,T,R}}, i, t) where {r,p,T,R}
# end
# function bsplinebasis₊₀(dP::BSplineDerivativeSpace{r,UniformBSplineSpace{p,T,R}}, i, t) where {r,p,T,R}
# end
# function bsplinebasis₋₀(dP::BSplineDerivativeSpace{r,UniformBSplineSpace{p,T,R}}, i, t) where {r,p,T,R}
# end
# function bsplinebasisall(dP::BSplineDerivativeSpace{r,UniformBSplineSpace{p,T,R}}, i, t) where {r,p,T,R}
# end
