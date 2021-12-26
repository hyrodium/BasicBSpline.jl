# Derivative

struct BSplineDerivativeSpace{r, T<:AbstractBSplineSpace}
    bsplinespace::T
    function BSplineDerivativeSpace{r}(P::T) where {r, T<:AbstractBSplineSpace}
        new{r,T}(P)
    end
end

bsplinespace(dP::BSplineDerivativeSpace) = dP.bsplinespace
knots(dP::BSplineDerivativeSpace) = knots(bsplinespace(dP))
degree(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = p - r
dim(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = dim(bsplinespace(dP)) - r

function _lower(dP::BSplineDerivativeSpace{r}) where r
    P = bsplinespace(dP)
    BSplineDerivativeSpace{r-1}(_lower(P))
end

@generated function bsplinebasis₊₀(P::BSplineDerivativeSpace{r,BSplineSpace{p,T}}, i::Integer, t::Real) where {r, p, T}
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:(T($(ks[i])≤t<$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one(T),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
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
            :(v = knots(P).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero(T))
    end
end

@generated function bsplinebasis₋₀(P::BSplineDerivativeSpace{r,BSplineSpace{p,T}}, i::Integer, t::Real) where {r, p, T}
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:(T($(ks[i])<t≤$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one(T),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
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
            :(v = knots(P).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero(T))
    end
end

@generated function bsplinebasis(P::BSplineDerivativeSpace{r,BSplineSpace{p,T}}, i::Integer, t::Real) where {r, p, T}
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:(T($(ks[i])≤t<$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    L_r(m,n) = Expr(:tuple, [:(_d(one(T),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
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
            :(v = knots(P).vector),
            :($k_l = $k_r),
            :($(B_l(p+1)) = $(A_r(p+1))),
            :($(Bs[end]) += $(T)(t == $(ks[end]) == v[end])),
            exs...,
            :(return $(prod(p-r+1:p))*B1)
        )
    else
        :(return zero(T))
    end
end

@doc raw"""
1st derivative of B-spline basis function.
Right-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
bsplinebasis′₊₀


@doc raw"""
1st derivative of B-spline basis function.
Left-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
bsplinebasis′₋₀

@doc raw"""
1st derivative of B-spline basis function.
Modified version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
bsplinebasis′

for suffix in ("", "₋₀", "₊₀")
    fname = Symbol(:bsplinebasis, suffix)
    fname′ = Symbol(:bsplinebasis, "′" ,suffix)
    @eval function $(fname′)(P::AbstractBSplineSpace, i::Integer, t::Real)
        return $(fname)(BSplineDerivativeSpace{1}(P), i, t)
    end
    @eval function $(fname′)(dP::BSplineDerivativeSpace{r}, i::Integer, t::Real) where r
        P = bsplinespace(dP)
        return $(fname)(BSplineDerivativeSpace{r+1}(P), i, t)
    end
end

@generated function bsplinebasisall(P::BSplineDerivativeSpace{r,BSplineSpace{p,T}}, i::Integer, t::Real) where {r, p, T}
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
            :(k = knots(P)),
            :($b = bsplinebasisall(_lower(P),i+1,t)),
            :($(Bs[1]) = $(K1s[1])*$(bs[1])),
            exs...,
            :($(Bs[p+1]) = $(K2s[p])*$(bs[p])),
            :(return $(B))
        )
    else
        Z = Expr(:tuple, [:(zero($T)) for i in 1:p+1]...)
        :(return $(Z))
    end
end

function bsplinebasisall(dP::BSplineDerivativeSpace{0,BSplineSpace{p,T}}, i::Integer, t::Real) where {p, T}
    P = bsplinespace(dP)
    bsplinebasisall(P,i,t)
end

# TODO: Add issubset(::BSplineDerivativeSpace, i, t)
# TODO: Add bsplineunity(::BSplineDerivativeSpace)
