# B-spline basis function

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
function bsplinebasis₊₀(P::BSplineSpace, t::Real)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] ≤ t < k[i+1] for i in 1:n]
    end
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₊₀(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
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
function bsplinebasis₋₀(P::BSplineSpace, t::Real)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] < t ≤ k[i+1] for i in 1:n]
    end
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₋₀(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
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
function bsplinebasis(P::BSplineSpace, t::Real)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t) for i in 1:n]
    end
    K = [ifelse(k[i+p] == k[i], 0, (t - k[i]) / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] + (1 - K[i+1]) * B[i+1] for i in 1:n]
end

"""
``i``-th B-spline basis function.
Right-sided limit version.
"""
function bsplinebasis₊₀(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    if p == 0
        return k[i] ≤ t < k[i+1]
    else
        return (
            ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₊₀(BSplineSpace(p-1, k), i, t) * (t - k[i]) / (k[i+p] - k[i]) : 0) +
            ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₊₀(BSplineSpace(p-1, k), i+1, t) * (k[i+p+1] - t) / (k[i+p+1] - k[i+1]) : 0)
        )
    end
end

"""
``i``-th B-spline basis function.
Left-sided limit version.
"""
function bsplinebasis₋₀(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    if p == 0
        return k[i] < t ≤ k[i+1]
    else
        return (
            ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₋₀(BSplineSpace(p-1, k), i, t) * (t - k[i]) / (k[i+p] - k[i]) : 0) +
            ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₋₀(BSplineSpace(p-1, k), i+1, t) * (k[i+p+1] - t) / (k[i+p+1] - k[i+1]) : 0)
        )
    end
end

"""
``i``-th B-spline basis function.
Modified version.
"""
function bsplinebasis(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    if p == 0
        return k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t)
    else
        return (
            ((k[i+p] - k[i] ≠ 0) ? bsplinebasis(BSplineSpace(p-1, k), i, t) * (t - k[i]) / (k[i+p] - k[i]) : 0) +
            ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis(BSplineSpace(p-1, k), i+1, t) * (k[i+p+1] - t) / (k[i+p+1] - k[i+1]) : 0)
        )
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
function bsplinebasis′₊₀(P::BSplineSpace, t)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₊₀(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
end

@doc raw"""
1st derivative of B-spline basis function.
Left-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
function bsplinebasis′₋₀(P::BSplineSpace, t)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis₋₀(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
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

function bsplinebasis′(P::BSplineSpace, t)::Vector{Float64}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p] == k[i], 0, p / (k[i+p] - k[i])) for i in 1:n+1]
    B = bsplinebasis(BSplineSpace(p-1, k), t)
    return [K[i] * B[i] - K[i+1] * B[i+1] for i in 1:n]
end

function bsplinebasis′₊₀(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₊₀(BSplineSpace(p-1, k), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₊₀(BSplineSpace(p-1, k), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end

function bsplinebasis′₋₀(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis₋₀(BSplineSpace(p-1, k), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis₋₀(BSplineSpace(p-1, k), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end

function bsplinebasis′(P::BSplineSpace, i::Integer, t::Real)::Float64
    p = P.degree
    k = P.knots

    return p * (
        ((k[i+p] - k[i] ≠ 0) ? bsplinebasis(BSplineSpace(p-1, k), i, t) / (k[i+p] - k[i]) : 0) -
        ((k[i+p+1] - k[i+1] ≠ 0) ? bsplinebasis(BSplineSpace(p-1, k), i+1, t) / (k[i+p+1] - k[i+1]) : 0)
    )
end

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```
"""
function bsplinesupport(P::AbstractBSplineSpace, i::Integer)
    p = degree(P)
    k = knots(P)
    return k[i]..k[i+p+1]
end

function bsplinesupport(P::AbstractBSplineSpace)
    p = degree(P)
    k = knots(P)
    return [k[i]..k[i+p+1] for i in 1:dim(P)]
end

# For compatibility
bsplinebasis₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinebasis′(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
bsplinesupport(i::Integer,P::AbstractBSplineSpace) = bsplinesupport(P,i)
