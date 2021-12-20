# B-Spline Basis Function

@inline _d(a::T,b::T) where T = ifelse(iszero(b), zero(T), T(a/b))

@doc raw"""
B-spline basis function.
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
```
"""
function bsplinebasis₊₀(P::BSplineSpace{p}, t::Real)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₊₀(lower(P), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
end
function bsplinebasis₊₀(P::BSplineSpace{0}, t::Real)::Vector{Float64}
    n = dim(P)
    k = knots(P)
    return [k[i] ≤ t < k[i+1] for i in 1:n]
end

@doc raw"""
B-spline basis function.
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
```
"""
function bsplinebasis₋₀(P::BSplineSpace{p}, t::Real)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] < t ≤ k[i+1] for i in 1:n]
    end
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₋₀(lower(P), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
end
function bsplinebasis₋₀(P::BSplineSpace{0}, t::Real)::Vector{Float64}
    n = dim(P)
    k = knots(P)
    return [k[i] < t ≤ k[i+1] for i in 1:n]
end

@doc raw"""
B-spline basis function.
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
```
"""
function bsplinebasis(P::BSplineSpace{p}, t::Real)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t) for i in 1:n]
    end
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis(lower(P), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
end
function bsplinebasis(P::BSplineSpace{0}, t::Real)::Vector{Float64}
    n = dim(P)
    k = knots(P)
    return [k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t) for i in 1:n]
end

"""
``i``-th B-spline basis function.
Right-sided limit version.
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

"""
``i``-th B-spline basis function.
Left-sided limit version.
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

"""
``i``-th B-spline basis function.
Modified version.
"""
function bsplinebasis(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
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
function bsplinebasis′₊₀(P::BSplineSpace{p}, t)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₊₀(lower(P), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
end
function bsplinebasis′₊₀(P::BSplineSpace{0}, t)::Vector{Float64}
    n = dim(P)
    return zeros(n)
end

@doc raw"""
1st derivative of B-spline basis function.
Left-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
function bsplinebasis′₋₀(P::BSplineSpace{p}, t)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₋₀(lower(P), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
end
function bsplinebasis′₋₀(P::BSplineSpace{0}, t)::Vector{Float64}
    n = dim(P)
    return zeros(n)
end

@doc raw"""
1st derivative of B-spline basis function.
Modified version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
bsplinebasis′

function bsplinebasis′(P::BSplineSpace{p}, t)::Vector{Float64} where p
    k = P.knots

    n = dim(P)
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis(lower(P), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
end
function bsplinebasis′(P::BSplineSpace{0}, t)::Vector{Float64}
    n = dim(P)
    return zeros(n)
end

function bsplinebasis′₊₀(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₊₀(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₊₀(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end

function bsplinebasis′₋₀(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₋₀(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₋₀(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end

function bsplinebasis′(P::BSplineSpace{p}, i::Integer, t::Real)::Float64 where p
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis(lower(P), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis(lower(P), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end


# For compatibility
bsplinebasis₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinesupport(i::Integer,P::AbstractBSplineSpace) = bsplinesupport(P,i)
