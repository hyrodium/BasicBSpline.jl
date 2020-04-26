# B-spline basis

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
function BSplineBasis₊₀(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] ≤ t < k[i+1] for i ∈ 1:n]
    end
    K = [ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis₊₀(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
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
function BSplineBasis₋₀(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] < t ≤ k[i+1] for i ∈ 1:n]
    end
    K = [ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis₋₀(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
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
function BSplineBasis(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return [k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t) for i ∈ 1:n]
    end
    K = [ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
end

"""
i-th B-spline basis function.
Modified version.
"""
function BSplineBasis(i::Int64, P::BSplineSpace, t)::Float64
    p = P.degree
    k = P.knots

    if p == 0
        return k[i]≤t<k[i+1] || (k[i]≠k[i+1]==k[end]==t)
    else
        return (((k[i+p]-k[i]≠0) ? BSplineBasis(i,𝒫(p-1,k),t)*(t-k[i])/(k[i+p]-k[i]) : 0)
        +((k[i+p+1]-k[i+1]≠0) ? BSplineBasis(i+1,𝒫(p-1,k),t)*(k[i+p+1]-t)/(k[i+p+1]-k[i+1]) : 0))
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
function BSplineBasis′₊₀(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis₊₀(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

@doc raw"""
1st derivative of B-spline basis function.
Left-sided limit version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
function BSplineBasis′₋₀(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis₋₀(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

@doc raw"""
1st derivative of B-spline basis function.
Modified version.
```math
\dot{B}_{(i,p,k)}(t)
=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
```
"""
function BSplineBasis′(P::BSplineSpace, t)::Array{Float64,1}
    p = P.degree
    k = P.knots

    n = dim(P)
    if p == 0
        return zeros(n)
    end
    K = [ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B = BSplineBasis(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

function BSplineBasis′(i::Int64, P::BSplineSpace, t)::Float64
    p = P.degree
    k = P.knots

    return p*(((k[i+p]-k[i]≠0) ? BSplineBasis(i,𝒫(p-1,k),t)/(k[i+p]-k[i]) : 0)
    -((k[i+p+1]-k[i+1]≠0) ? BSplineBasis(i+1,𝒫(p-1,k),t)/(k[i+p+1]-k[i+1]) : 0))
end

@doc raw"""
Return support of i-th B-spline basis function.
"""
function BSplineSupport(i::Int64, P::BSplineSpace)::ClosedInterval
    p = P.degree
    k = P.knots
    return k[i]..k[i+p+1]
end
