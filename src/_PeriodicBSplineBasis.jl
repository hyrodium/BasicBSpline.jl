######################################
#= Periodic B-Spline Basis Function =#
######################################

@doc raw"""
Return the support of the ``i``-th periodic B-spline basis function, represented on the
cyclically-extended knot axis:
```math
\operatorname{supp}(B_{(i,p,k)}) = [\tilde{k}_i, \tilde{k}_{i+p+1}]
```
where ``\tilde{k}_j = k_{((j-1) \bmod n) + 1} + L \lfloor (j-1)/n \rfloor``.

The returned interval may extend past one fundamental period; taking modulo the period
gives the actual (possibly wrapped) support on the circle.
"""
function bsplinesupport_R(P::PeriodicBSplineSpace{p}, i::Integer) where p
    return _periodic_knot(P, i)..(_periodic_knot(P, i+p+1))
end

@doc raw"""
Periodic B-spline basis functions have no boundary on the circle; `bsplinesupport_I`
coincides with `bsplinesupport_R`.
"""
bsplinesupport_I(P::PeriodicBSplineSpace, i::Integer) = bsplinesupport_R(P, i)

"""
Return the fundamental interval index containing `t` (after reduction modulo the period).
Result is in `1:n` where `n = dim(P)`.
"""
function intervalindex(P::PeriodicBSplineSpace, t::Real)
    k = knotvector(P)
    L = period(P)
    t0 = k[1]
    t_red = mod(t - t0, L) + t0
    v = _vec(k)
    j = searchsortedlast(v, t_red)
    n = length(k) - 1
    return clamp(j, 1, n)
end

@doc raw"""
``i``-th periodic B-spline basis function.
Right-sided limit version.

Computed by the Cox–de Boor recurrence applied to the cyclically-extended knot sequence,
with `t` shifted into the support interval ``[\tilde{k}_i, \tilde{k}_i + L)``.
"""
@generated function bsplinebasis₊₀(P::PeriodicBSplineSpace{p,T}, i::Integer, t::S) where {p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    ks = [Symbol(:k,j) for j in 1:p+2]
    Ks = [Symbol(:K,j) for j in 1:p+1]
    Bs = [Symbol(:B,j) for j in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, [:(_periodic_knot(P, i+$(j-1))) for j in 1:p+2]...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:($U($(ks[j])≤ts<$(ks[j+1]))) for j in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(ts-$(ks[j]),$(ks[j+m])-$(ks[j]))) for j in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[j])*$(Bs[j])+(1-$(Ks[j+1]))*$(Bs[j+1])) for j in 1:n]...)
    exs = Expr[]
    for j in 1:p
        push!(exs, :($(K_l(p+2-j)) = $(K_r(j,p+2-j))))
        push!(exs, :($(B_l(p+1-j)) = $(B_r(p+1-j))))
    end
    return Expr(:block,
        :(L = period(P)),
        :($k_l = $k_r),
        :(ts = mod(t - $(ks[1]), L) + $(ks[1])),
        :($(B_l(p+1)) = $(A_r(p+1))),
        exs...,
        :(return B1)
    )
end

@doc raw"""
``i``-th periodic B-spline basis function.
Left-sided limit version.
"""
@generated function bsplinebasis₋₀(P::PeriodicBSplineSpace{p,T}, i::Integer, t::S) where {p, T, S<:Real}
    U = StaticArrays.arithmetic_closure(promote_type(T,S))
    ks = [Symbol(:k,j) for j in 1:p+2]
    Ks = [Symbol(:K,j) for j in 1:p+1]
    Bs = [Symbol(:B,j) for j in 1:p+1]
    k_l = Expr(:tuple, ks...)
    k_r = Expr(:tuple, [:(_periodic_knot(P, i+$(j-1))) for j in 1:p+2]...)
    K_l(n) = Expr(:tuple, Ks[1:n]...)
    B_l(n) = Expr(:tuple, Bs[1:n]...)
    A_r(n) = Expr(:tuple, [:($U($(ks[j])<ts≤$(ks[j+1]))) for j in 1:n]...)
    K_r(m,n) = Expr(:tuple, [:(_d(ts-$(ks[j]),$(ks[j+m])-$(ks[j]))) for j in 1:n]...)
    B_r(n) = Expr(:tuple, [:($(Ks[j])*$(Bs[j])+(1-$(Ks[j+1]))*$(Bs[j+1])) for j in 1:n]...)
    exs = Expr[]
    for j in 1:p
        push!(exs, :($(K_l(p+2-j)) = $(K_r(j,p+2-j))))
        push!(exs, :($(B_l(p+1-j)) = $(B_r(p+1-j))))
    end
    return Expr(:block,
        :(L = period(P)),
        :($k_l = $k_r),
        :(ts = mod(t - $(ks[1]), L) + $(ks[1])),
        :($(B_l(p+1)) = $(A_r(p+1))),
        exs...,
        :(return B1)
    )
end

@doc raw"""
``i``-th periodic B-spline basis function. On the circle the left and right limits
coincide, so this is identical to [`bsplinebasis₊₀`](@ref).
"""
@inline bsplinebasis(P::PeriodicBSplineSpace, i::Integer, t::Real) = bsplinebasis₊₀(P, i, t)
