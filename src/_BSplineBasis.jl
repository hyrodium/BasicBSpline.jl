#############################
#= B-Spline Basis Function =#
#############################

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

# bsplinebasis with UniformKnotVector
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

################################################
#= B-Spline Basis Functions at specific point =#
################################################

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

@generated function bsplinebasisall(P::BSplineSpace{p,T}, i::Integer, t::Real) where {p,T}
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

## bsplinebasisall with UniformKnotVector
@inline function bsplinebasisall(P::BSplineSpace{p,T,<:UniformKnotVector{T,R}},i::Integer,t::S) where {p, T, R<:AbstractUnitRange, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (0 ≤ i ≤ length(k)-2p) || throw(DomainError(i, "index of interval is out of range."))
    return uniform_bsplinebasisall_kernel(Val{p}(),U(t-i+1-p-k[1]))
end
@inline function bsplinebasisall(P::BSplineSpace{p,T,<:UniformKnotVector{T}},i::Integer,t::S) where {p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    k = knotvector(P)
    @boundscheck (0 ≤ i ≤ length(k)-2p) || throw(DomainError(i, "index of interval is out of range."))
    a = @inbounds k.vector[i+p]
    b = @inbounds k.vector[i+p+1]
    return uniform_bsplinebasisall_kernel(Val{p}(),U((t-a)/(b-a)))
end


#########################################################
#= Kernel funcitons of uniform B-spline basis function =#
#########################################################

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
