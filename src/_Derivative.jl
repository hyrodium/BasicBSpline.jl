# Derivative

struct BSplineDerivativeSpace{r, T<:AbstractBSplineSpace}
    bsplinespace::T
    function BSplineDerivativeSpace{r}(P::T) where {r, T<:AbstractBSplineSpace}
        new{r,T}(P)
    end
end

bsplinespace(P::BSplineDerivativeSpace) = P.bsplinespace

knots(P::BSplineDerivativeSpace) = knots(bsplinespace(P))

degree(P::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = p - r

dim(P::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}) where {r,p} = dim(bsplinespace(P)) - r

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

# TODO: Add bsplinebasisall(::BSplineDerivativeSpace, i, t)
# TODO: Add issubset(::BSplineDerivativeSpace, i, t)
# TODO: Add bsplineunity(::BSplineDerivativeSpace)
