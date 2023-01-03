# B-Spline Basis Function

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
```

# Examples
```jldoctest
julia> P = BSplineSpace{0}(KnotVector(1:6))
BSplineSpace{0, Int64, KnotVector{Int64}}(KnotVector([1, 2, 3, 4, 5, 6]))

julia> bsplinebasis₊₀.(P,1:5,(1:6)')
5×6 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
```
"""
@generated function bsplinebasis₊₀(P::BSplineSpace{p,T}, i::Integer, t::S) where {p, T, S<:Real}
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
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    exs = Expr[]
    for i in 1:p
        push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
        push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
    end
    Expr(:block,
        :(v = knotvector(P).vector),
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
```

# Examples
```jldoctest
julia> P = BSplineSpace{0}(KnotVector(1:6))
BSplineSpace{0, Int64, KnotVector{Int64}}(KnotVector([1, 2, 3, 4, 5, 6]))

julia> bsplinebasis₋₀.(P,1:5,(1:6)')
5×6 Matrix{Float64}:
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0
```
"""
@generated function bsplinebasis₋₀(P::BSplineSpace{p,T}, i::Integer, t::S) where {p, T, S<:Real}
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
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    exs = Expr[]
    for i in 1:p
        push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
        push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
    end
    Expr(:block,
        :(v = knotvector(P).vector),
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
    &1\quad (k_{i} \le t<k_{i+1})\\
    &1\quad (k_{i} < t = k_{i+1}=k_{l})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

# Examples
```jldoctest
julia> P = BSplineSpace{0}(KnotVector(1:6))
BSplineSpace{0, Int64, KnotVector{Int64}}(KnotVector([1, 2, 3, 4, 5, 6]))

julia> bsplinebasis.(P,1:5,(1:6)')
5×6 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  1.0
```
"""
@generated function bsplinebasis(P::BSplineSpace{p,T}, i::Integer, t::S) where {p, T, S<:Real}
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
    B_r(n) = Expr(:tuple, [:($(Ks[i])*$(Bs[i])+(1-$(Ks[i+1]))*$(Bs[i+1])) for i in 1:n]...)
    exs = Expr[]
    for i in 1:p
        push!(exs, :($(K_l(p+2-i)) = $(K_r(i,p+2-i))))
        push!(exs, :($(B_l(p+1-i)) = $(B_r(p+1-i))))
    end
    Expr(:block,
        :(v = knotvector(P).vector),
        :($k_l = $k_r),
        :($(B_l(p+1)) = $(A_r(p+1))),
        exs...,
        :(return B1)
    )
end

"""
B-spline basis functions at point `t` on `i`-th interval.

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

julia> p = 2
2

julia> P = BSplineSpace{p}(k)
BSplineSpace{2, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> t = 5.7
5.7

julia> i = intervalindex(P,t)
2

julia> bsplinebasisall(P,i,t)
3-element SVector{3, Float64} with indices SOneTo(3):
 0.3847272727272727
 0.6107012987012989
 0.00457142857142858

julia> bsplinebasis.(P,i:i+p,t)
3-element Vector{Float64}:
 0.38472727272727264
 0.6107012987012989
 0.00457142857142858
```
"""
bsplinebasisall

@inline function bsplinebasisall(P::BSplineSpace{0,T},i::Integer,t::S) where {T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    SVector(one(U),)
end

@inline function bsplinebasisall(P::BSplineSpace{1}, i::Integer, t::Real)
    k = knotvector(P)
    B1 = (k[i+2]-t)/(k[i+2]-k[i+1])
    B2 = (t-k[i+1])/(k[i+2]-k[i+1])
    return SVector(B1, B2)
end

@generated function bsplinebasisall(P::BSplineSpace{p}, i::Integer, t::Real) where p
    bs = [Symbol(:b,i) for i in 1:p]
    Bs = [Symbol(:B,i) for i in 1:p+1]
    K1s = [:((k[i+$(p+j)]-t)/(k[i+$(p+j)]-k[i+$(j)])) for j in 1:p]
    K2s = [:((t-k[i+$(j)])/(k[i+$(p+j)]-k[i+$(j)])) for j in 1:p]
    b = Expr(:tuple, bs...)
    B = Expr(:tuple, Bs...)
    exs = [:($(Bs[j+1]) = ($(K1s[j+1])*$(bs[j+1]) + $(K2s[j])*$(bs[j]))) for j in 1:p-1]
    Expr(:block,
        :($(Expr(:meta, :inline))),
        :(k = knotvector(P)),
        :($b = bsplinebasisall(_lower(P),i+1,t)),
        :($(Bs[1]) = $(K1s[1])*$(bs[1])),
        exs...,
        :($(Bs[p+1]) = $(K2s[p])*$(bs[p])),
        :(return SVector($(B)))
    )
end
