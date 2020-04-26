module BasicBSpline

using IntervalSets

export Knots, BSplineSpace, 𝒫, dim
export BSplineBasis₊₀, BSplineBasis₋₀, BSplineBasis
export BSplineBasis′₊₀, BSplineBasis′₋₀, BSplineBasis′
export BSplineSupport, BSplineCoefficient
export BSplineManifold, Refinement, Mapping, BSplineSvg

# Knots
@doc raw"""
Construct knot vector from given array.
```math
k=(k_1,\dots,k_l)
```
"""
struct Knots
    vector :: Array{Float64,1}
    function Knots(vector::AbstractArray{T,1} where T<:Real)
        return new(sort(convert(Array{Float64,1},vector)))
    end
    function Knots(vector::Array{Any,1})
        if isempty(vector)
            return new(Float64[])
        else
            error("The elements of given vector must be real number.")
        end
    end
end

Base.zero(::Type{Knots}) = Knots([])
Base. ==(k₁::Knots, k₂::Knots) = (k₁.vector==k₂.vector)
Base.:+(k₁::Knots, k₂::Knots) = Knots(sort([k₁.vector...,k₂.vector...]))
Base.:*(p₊::Int, k::Knots) = (
        if p₊ == 0
            zero(Knots)
        elseif p₊ > 0
            sum(k for _ ∈ 1:p₊)
        else
            error("Polynominal degree p₊ must be non-negative.")
        end
    )

Base.in(r::Real, k::Knots) = in(r,k.vector)
Base.getindex(k::Knots, i::Int) = k.vector[i]
Base.getindex(k::Knots, v::AbstractArray{Int64,1}) = Knots(k.vector[v])
Base.length(k::Knots) = length(k.vector)
♯(k::Knots) = length(k::Knots)
Base.firstindex(k) = 1
Base.lastindex(k) = length(k)
Base.unique(k::Knots) = Knots(unique(k.vector))

function Base.:⊆(k::Knots, k′::Knots)
    K′ = copy(k′.vector)
    for kᵢ ∈ k.vector
        i = findfirst(x -> x == kᵢ,K′)
        if i isa Nothing
            return false
        end
        deleteat!(K′,i)
    end
    return true
end


# B-Spline Space
@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```
"""
struct BSplineSpace
    degree::Int
    knots::Knots
    function BSplineSpace(degree::Int, knots::Knots)
        if degree < 0
            error("degree of polynominal must be non-negative")
        end
        new(degree,knots)
    end
end

@doc raw"""
Same as BSplineSpace.
```math
\mathcal{P}[p,k]
```
"""
const 𝒫 = BSplineSpace

@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\sharp k - p -1
```
"""
function dim(bsplinespace::BSplineSpace)
    p=bsplinespace.degree
    k=bsplinespace.knots
    return ♯(k)-p-1
end

"""
Check inclusive relationship between B-spline spaces.
"""
function Base.:⊆(P::BSplineSpace, P′::BSplineSpace)
    p = P.degree
    k = P.knots
    p′ = P′.degree
    k′ = P′.knots
    p₊ = p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function Base.iszero(P::BSplineSpace)
    p = P.degree
    k = P.knots
    n = dim(P)
    return [k[i] == k[i+p+1] for i ∈ 1:n]
end

# B-Spline functions
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

function BSplineCoefficient(P::BSplineSpace, P′::BSplineSpace)::Array{Float64,2}
    p = P.degree
    k = P.knots
    p′ = P′.degree
    k′ = P′.knots
    p₊ = p′-p
    if P ⊈ P′
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end

    if p == 0
        n=length(k)-1
        n′=length(k′)-p₊-1
        A⁰=Float64[BSplineSupport(j,𝒫(p₊,k′)) ⊆ BSplineSupport(i,𝒫(0,k)) for i ∈ 1:n, j ∈ 1:n′]
        A⁰[:,findall(iszero(P′))].=NaN
        return A⁰
    end

    Aᵖ⁻¹=BSplineCoefficient(𝒫(p-1, k), 𝒫(p′-1, k′))
    n = dim(P)
    n′=dim(P′)
    Z = iszero(𝒫(p′-1,k′))
    W = findall(Z)
    K′ = [k′[i+p′]-k′[i] for i ∈ 1:n′+1]
    K = [ifelse(k[i+p]≠k[i], 1/(k[i+p]-k[i]), 0.0) for i ∈ 1:n+1]
    Δ = (p/p′)*[K′[j]*(K[i]*Aᵖ⁻¹[i,j]-K[i+1]*Aᵖ⁻¹[i+1,j]) for i ∈ 1:n, j ∈ 1:n′+1]
    Aᵖ = zeros(n,n′)
    Aᵖ[:,1] = Δ[:,1]
    Aᵖ[:,n′] = -Δ[:,n′+1]

    if length(W) == 0
        Q = [1:n′]
    else
        Q = [1:W[1]-1,[W[i]:W[i+1]-1 for i ∈ 1:length(W)-1]...,W[end]:n′]
    end
    l = length(Q)
    L = length.(Q)
    Ãᵖ = [Aᵖ[:,q] for q ∈ Q]

    for ȷ ∈ 2:l-1
        if L[ȷ] == 1
            Ãᵖ[ȷ] .= NaN
        end
    end
    for ȷ ∈ 1:l-1
        if L[ȷ] ≥ 2
            t = k′[W[ȷ]]
            Ãᵖ[ȷ][:,end] = BSplineBasis₋₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 2:l
        if L[ȷ] ≥ 2
            t = k′[W[ȷ-1]+p]
            Ãᵖ[ȷ][:,1] = BSplineBasis₊₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 1:l
        if L[ȷ] ≥ 3
            r = Q[ȷ]
            A₊ = copy(Ãᵖ[ȷ])
            A₋ = copy(Ãᵖ[ȷ])
            for j ∈ 1:L[ȷ]-2
                A₊[:,j+1] = A₊[:,j]+Δ[:,j+r[1]]
                A₋[:,L[ȷ]-j] = A₋[:,L[ȷ]-j+1]-Δ[:,L[ȷ]-j+r[1]]
            end
            Ãᵖ[ȷ] = (A₊+A₋)/2
        end
    end
    Aᵖ = hcat(Ãᵖ...)
    return Aᵖ .* Float64[BSplineSupport(j,P′) ⊆ BSplineSupport(i,P) for i ∈ 1:n, j ∈ 1:n′]
end

function ⊗(X::Array{Float64},Y::Array{Float64})::Array{Float64}
    m = size(X)
    n = size(Y)
    reshape(reshape(X,length(X)) * reshape(Y,length(Y))', m..., n...)
end

function tensorprod(X::Array{T,1}) where T <: Array{Float64}
    n=length(X)
    # X[1] ⊗ … ⊗ X[n]
    @inbounds Y = X[1]
    for i ∈ 2:n
        @inbounds Y = Y ⊗ X[i]
    end
    return Y
end

struct BSplineManifold
    bsplinespaces::Array{BSplineSpace,1}
    controlpoints::Array{Float64}
    function BSplineManifold(bsplinespaces::Array{BSplineSpace,1}, controlpoints::Array{Float64})
        if collect(size(controlpoints)[1:end-1]) ≠ dim.(bsplinespaces)
            error("dimension does not match")
        else
            new(bsplinespaces, controlpoints)
        end
    end
end

@doc raw"""
Refinement of B-spline manifold.
"""
function Refinement(M::BSplineManifold, Ps′::Array{BSplineSpace,1})
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d̂ = size(𝒂)[end]
    n = dim.(Ps)
    n′ = dim.(Ps′)
    if prod(Ps .⊆ Ps′)
        A = BSplineCoefficient.(Ps,Ps′)
        𝒂′ = [sum(A[1][I₁,J₁]*A[2][I₂,J₂]*𝒂[I₁,I₂,i] for I₁ ∈ 1:n[1], I₂ ∈ 1:n[2]) for J₁ ∈ 1:n′[1], J₂ ∈ 1:n′[2], i ∈ 1:d̂]
        return BSplineManifold(Ps′, 𝒂′)
    else
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end
end

@doc raw"""
Refinement of B-spline manifold.
"""
function Refinement(M::BSplineManifold; p₊::Union{Nothing,Array{Int,1}}=nothing, k₊::Union{Nothing,Array{Knots,1}}=nothing)
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d = length(Ps)
    d̂ = size(𝒂)[end]
    n = dim.(Ps)
    if p₊ == nothing
        p₊=zeros(Int,d)
    elseif length(Ps) ≠ length(p₊)
        error("dimension does not match")
    end
    if k₊ == nothing
        k₊=zeros(Knots,d)
    elseif length(Ps) ≠ length(k₊)
        error("dimension does not match")
    end

    Ps′ = BSplineSpace[]
    for i ∈ 1:length(Ps)
        P = Ps[i]
        p = P.degree
        k = P.knots
        push!(Ps′,𝒫(p+p₊[i], k+p₊[i]*unique(k)+k₊[i]))
    end

    return Refinement(M, Ps′)
end

# function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
#     if length(𝒫s)==length(t)==1
#         return BSplineBasis(𝒫s[1],t[1])
#     elseif length(𝒫s)==length(t)==2
#         return BSplineBasis(𝒫s[1],t[1])*BSplineBasis(𝒫s[2],t[2])'
#     else
#         error("dimension does not match")
#     end
# end

@doc raw"""
Multi-dimentional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
    d = length(t)
    Bs = [BSplineBasis(𝒫s[i],t[i]) for i ∈ 1:d]
    return tensorprod(Bs)
end

# function Mapping(M::BSplineManifold, t::Array{Float64,1})
#     𝒫s = M.bsplinespaces
#     𝒂 = M.controlpoints
#     d=length(𝒫s)
#     d̂=size(𝒂)[end]
#     return [sum(BSplineBasis(𝒫s,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
# end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function Mapping(M::BSplineManifold, t::Array{Float64,1})
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d = length(Ps)
    d̂ = size(𝒂)[end]
    return [sum(BSplineBasis(Ps,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
end

end # module
