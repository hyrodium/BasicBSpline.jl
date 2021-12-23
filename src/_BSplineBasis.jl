# B-Spline Basis Function

@inline _d(a::T,b::T) where T = ifelse(iszero(b), zero(T), T(a/b))

@doc raw"""
``i``-th B-spline basis function.
Right-sided limit version.
```math
\begin{aligned}
{B}_{(i,p,k)}(t)
&=
\frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
+\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}\le t< k_{i+1})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
"""
@generated function bsplinebasis₊₀(P::BSplineSpace{p,T}, i::Integer, t::Real) where {p, T}
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:(T($(ks[i])≤t<$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    exs = Expr[]
    for i in 1:p
        push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
        push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
    end
    Expr(:block,
        :(v = knots(P).vector),
        :($k_l = $k_r),
        :($(B_l(p+1)) = $(A_r(p+1))),
        exs...,
        :(return B1)
    )
end

@doc raw"""
``i``-th B-spline basis function.
Left-sided limit version.
```math
\begin{aligned}
{B}_{(i,p,k)}(t)
&=
\frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
+\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}< t\le k_{i+1})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
"""
@generated function bsplinebasis₋₀(P::BSplineSpace{p,T}, i::Integer, t::Real) where {p, T}
    ks = [Symbol(:k,i) for i in 1:p+2]
    Ks = [Symbol(:K,i) for i in 1:p+1]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, :(v[i]), (:(v[i+$j]) for j in 1:p+1)...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:(T($(ks[i])<t≤$(ks[i+1]))) for i in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(t-$(ks[i]),$(ks[i+m])-$(ks[i]))) for i in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    exs = Expr[]
    for i in 1:p
        push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
        push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
    end
    Expr(:block,
        :(v = knots(P).vector),
        :($k_l = $k_r),
        :($(B_l(p+1)) = $(A_r(p+1))),
        exs...,
        :(return B1)
    )
end

@doc raw"""
``i``-th B-spline basis function.
Modified version.
```math
\begin{aligned}
{B}_{(i,p,k)}(t)
&=
\frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
+\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}\le t<k_{i+1}<k_{l})\\
    &1\quad (k_{i}\le t\le k_{i+1}=k_{l})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
"""
function bsplinebasis(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    # TODO: use @generated macro
    k = P.knots

    return (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis(lower(P), i, t) * (t - k[i]) / (k[i+p] - k[i]) : 0) +
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis(lower(P), i+1, t) * (k[i+p+1] - t) / (k[i+p+1] - k[i+1]) : 0)
    )
end
function bsplinebasis(P::BSplineSpace{0}, i::Integer, t::Real)::Float64
    k = knots(P)
    return k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t)
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

function bsplinebasis′₊₀(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots
    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₊₀(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₊₀(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end
function bsplinebasis′₊₀(P::BSplineSpace{0}, i::Integer, t::Real)::Float64
    return 0.0
end

function bsplinebasis′₋₀(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots
    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₋₀(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₋₀(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end
function bsplinebasis′₋₀(P::BSplineSpace{0}, i::Integer, t::Real)::Float64
    return 0.0
end

function bsplinebasis′(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots
    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end
function bsplinebasis′(P::BSplineSpace{0}, i::Integer, t::Real)::Float64
    return 0.0
end


"""
TODO: Add docstring
"""
bsplinebasisall

# TODO: faster implementation
# TODO: use @generated macro
function bsplinebasisall(P::BSplineSpace{0,T},i::Integer,t::T) where T
    (one(T),)
end

function bsplinebasisall(P::BSplineSpace{1,T}, i::Integer, t::T) where T
    k = knots(P)
    B1 = (k[i+2]-t)/(k[i+2]-k[i+1])
    B2 = (t-k[i+1])/(k[i+2]-k[i+1])
    return (B1, B2)
end

function bsplinebasisall(P::BSplineSpace{2,T}, i::Integer, t::T) where T
    k = knots(P)
    B = bsplinebasisall(lower(P),i+1,t)

    B1 = (k[i+3]-t)/(k[i+3]-k[i+1]) * B[1]
    B2 = (t-k[i+1])/(k[i+3]-k[i+1]) * B[1] +
         (k[i+4]-t)/(k[i+4]-k[i+2]) * B[2]
    B3 = (t-k[i+2])/(k[i+4]-k[i+2]) * B[2]
    return B1, B2, B3
end

function bsplinebasisall(P::BSplineSpace{3,T}, i::Integer, t::T) where T
    k = knots(P)
    B = bsplinebasisall(lower(P),i+1,t)

    B1 = (k[i+4]-t)/(k[i+4]-k[i+1]) * B[1]
    B2 = (t-k[i+1])/(k[i+4]-k[i+1]) * B[1] +
         (k[i+5]-t)/(k[i+5]-k[i+2]) * B[2]
    B3 = (t-k[i+2])/(k[i+5]-k[i+2]) * B[2] +
         (k[i+6]-t)/(k[i+6]-k[i+3]) * B[3]
    B4 = (t-k[i+3])/(k[i+6]-k[i+3]) * B[3]
    return B1, B2, B3, B4
end

function bsplinebasisall(P::BSplineSpace{4,T}, i::Integer, t::T) where T
    k = knots(P)
    B = bsplinebasisall(lower(P),i+1,t)

    B1 = (k[i+5]-t)/(k[i+5]-k[i+1]) * B[1]
    B2 = (t-k[i+1])/(k[i+5]-k[i+1]) * B[1] +
         (k[i+6]-t)/(k[i+6]-k[i+2]) * B[2]
    B3 = (t-k[i+2])/(k[i+6]-k[i+2]) * B[2] +
         (k[i+7]-t)/(k[i+7]-k[i+3]) * B[3]
    B4 = (t-k[i+3])/(k[i+7]-k[i+3]) * B[3] +
         (k[i+8]-t)/(k[i+8]-k[i+4]) * B[4]
    B5 = (t-k[i+4])/(k[i+8]-k[i+4]) * B[4]
    return B1, B2, B3, B4, B5
end

function bsplinebasisall(P::BSplineSpace{5,T}, i::Integer, t::T) where T
    k = knots(P)
    B = bsplinebasisall(lower(P),i+1,t)

    B1 = (k[i+6]-t)/(k[i+6]-k[i+1]) * B[1]
    B2 = (t-k[i+1])/(k[i+6]-k[i+1]) * B[1] +
         (k[i+7]-t)/(k[i+7]-k[i+2]) * B[2]
    B3 = (t-k[i+2])/(k[i+7]-k[i+2]) * B[2] +
         (k[i+8]-t)/(k[i+8]-k[i+3]) * B[3]
    B4 = (t-k[i+3])/(k[i+8]-k[i+3]) * B[3] +
         (k[i+9]-t)/(k[i+9]-k[i+4]) * B[4]
    B5 = (t-k[i+4])/(k[i+9]-k[i+4]) * B[4] +
         (k[i+10]-t)/(k[i+10]-k[i+5]) * B[5]
    B6 = (t-k[i+5])/(k[i+10]-k[i+5]) * B[5]
    return B1, B2, B3, B4, B5, B6
end
